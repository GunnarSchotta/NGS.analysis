#!/usr/bin/env python3
"""
NGS analysis
"""

__auhor__ = ["Gunnar Schotta"]
__email__ = "gunnar.schotta@bmc.med.lmu.de"
__version__ = "0.9.5"


from argparse import ArgumentParser
import os
import re
import sys
import tempfile
import subprocess
import yaml
import pypiper
from pypiper import build_command

PROTOCOLS = ["CHIP", "CT", "CR", "RNA", "ATAC"]
PIPELINEMODE = ["genes", "repeats"]

"""
Main pipeline process.
"""
parser = ArgumentParser(description='NGS analysis version ' + __version__)
parser = pypiper.add_pypiper_args(parser, groups=
    ['pypiper', 'looper', 'ngs'],
    required=["input", "genome", "sample-name", "output-parent"])
parser.add_argument("--protocol", dest="protocol", type=str,
                        default=None, choices=PROTOCOLS,
                        help="NGS processing protocol.")
parser.add_argument("--pipeline-mode", dest="pipeline_mode", type=str,
                        default=None, choices=PIPELINEMODE,
                        help="Pipeline run mode.")
parser.add_argument("--STAR_RNA_index", dest="STAR_RNA_index", type=str,
                        default=None,
                        help="STAR RNA index.")
parser.add_argument("--STAR_genome_index", dest="STAR_genome_index", type=str,
                        default=None,
                        help="STAR genome index.")
parser.add_argument("--Bowtie2_index", dest="Bowtie2_index", type=str,
                        default=None,
                        help="Bowtie2 genome index.")
parser.add_argument("--refgene_tss", dest="refgene_tss", type=str,
                        default=None,
                        help="Refgene TSS Bed file")
parser.add_argument("--genome_index", dest="genome_index", type=str,
                        default=None,
                        help="Genome Fasta index file (.fai)")
parser.add_argument("--repeats_SAF", dest="repeats_SAF", type=str,
                        default=None,
                        help="Repeats SAF file for featureCounts")
parser.add_argument("--repeats_SAFid", dest="repeats_SAFid", type=str,
                        default=None,
                        help="Repeats SAF file for featureCounts, individual Repeats")

args = parser.parse_args()
args.paired_end = args.single_or_paired.lower() == "paired"

if not args.input:
    parser.print_help()
    raise SystemExit

# Initialize, creating global PipelineManager and NGSTk instance for
# access in ancillary functions outside of main().
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

global pm
pm = pypiper.PipelineManager(
    name="NGS analysis", outfolder=outfolder, args=args, version=__version__)

global ngstk
ngstk = pypiper.NGSTk(pm=pm)

# Convenience alias
tools = pm.config.tools
param = pm.config.parameters
res = pm.config.resources

def tool_path(tool_name):
    """
    Return the path to a tool used by this pipeline.

    :param str tool_name: name of the tool (e.g., a script file_name)
    :return str: real, absolute path to tool (expansion and symlink resolution)
    """

    return os.path.join(os.path.dirname(__file__), tool_name)


############################################################################
#          Check that the input file(s) exist before continuing            #
############################################################################
if os.path.isfile(args.input[0]) and os.stat(args.input[0]).st_size > 0:
    print("Local input file: " + args.input[0])
elif os.path.isfile(args.input[0]) and os.stat(args.input[0]).st_size == 0:
    # The read1 file exists but is empty
    err_msg = "File exists but is empty: {}"
    pm.fail_pipeline(IOError(err_msg.format(args.input[0])))
else:
    # The read1 file does not exist
    err_msg = "Could not find: {}"
    pm.fail_pipeline(IOError(err_msg.format(args.input[0])))

if args.input2:
    if (os.path.isfile(args.input2[0]) and
            os.stat(args.input2[0]).st_size > 0):
        print("Local input file: " + args.input2[0])
    elif (os.path.isfile(args.input2[0]) and
            os.stat(args.input2[0]).st_size == 0):
        # The read1 file exists but is empty
        err_msg = "File exists but is empty: {}"
        pm.fail_pipeline(IOError(err_msg.format(args.input2[0])))
    else:
        # The read1 file does not exist
        err_msg = "Could not find: {}"
        pm.fail_pipeline(IOError(err_msg.format(args.input2[0])))

container = None  # legacy

############################################################################
#                      Grab and prepare input files                        #
############################################################################
pm.report_result(
    "File_mb",
    round(ngstk.get_file_size(
          [x for x in [args.input, args.input2] if x is not None])), 2)
pm.report_result("Read_type", args.single_or_paired)
pm.report_result("Genome", args.genome_assembly)

############################################################################
# 				Analysis pipeline 			   #
############################################################################

# Each (major) step should have its own subfolder
raw_folder = os.path.join(outfolder, "raw")
fastq_folder = os.path.join(outfolder, "fastq")

pm.timestamp("### Merge/link and fastq conversion: ")
# This command will merge multiple inputs so you can use multiple
# sequencing lanes in a single pipeline run.
local_input_files = ngstk.merge_or_link(
    [args.input, args.input2], raw_folder, args.sample_name)
# flatten nested list
if any(isinstance(i, list) for i in local_input_files):
    local_input_files = [i for e in local_input_files for i in e]
# maintain order and remove duplicate entries
local_input_files = list(dict.fromkeys(local_input_files))

cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(
    local_input_files, args.sample_name, args.paired_end, fastq_folder,
    zipmode=True)

# flatten nested list
if any(isinstance(i, list) for i in unaligned_fastq):
    unaligned_fastq = [i for e in unaligned_fastq for i in e]
# maintain order and remove duplicate entries
if any(isinstance(i, dict) for i in local_input_files):
    unaligned_fastq = list(dict.fromkeys(unaligned_fastq))

pm.run(cmd, unaligned_fastq,
       follow=ngstk.check_fastq(
        local_input_files, unaligned_fastq, args.paired_end))
pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)

if args.paired_end:
    untrimmed_fastq1 = unaligned_fastq[0]
    untrimmed_fastq2 = unaligned_fastq[1]
else:
    untrimmed_fastq1 = unaligned_fastq
    untrimmed_fastq2 = None

# Also run a fastqc (if installed/requested)
fastqc_folder = os.path.join(outfolder, "fastqc")
fastqc_report = os.path.join(fastqc_folder,
    args.sample_name + "_R1_fastqc.html")
fastqc_report_R2 = os.path.join(fastqc_folder,
    args.sample_name + "_R2_fastqc.html")
if ngstk.check_command(tools.fastqc):
        ngstk.make_dir(fastqc_folder)
if fastqc_folder and os.path.isabs(fastqc_folder):
    ngstk.make_sure_path_exists(fastqc_folder)
cmd = (tools.fastqc + " --noextract --outdir " +
       fastqc_folder + " " + untrimmed_fastq1)
pm.run(cmd, fastqc_report, nofail=False)
pm.report_object("FastQC report r1", fastqc_report)

if args.paired_end and untrimmed_fastq2:
    cmd = (tools.fastqc + " --noextract --outdir " +
           fastqc_folder + " " + untrimmed_fastq2)
    pm.run(cmd, fastqc_report_R2, nofail=False)
    pm.report_object("FastQC report r2", fastqc_report_R2)


############################################################################
#                     Adapter trimming  		                   #
############################################################################

pm.timestamp("### Adapter trimming: ")

#Adapters: Nextera for ATAC and C&T; Illumina for rest
if args.protocol in ["ATAC", "CT"]:
    adapters = res.adapters_nextera
else:
    adapters = res.adapters_illumina

# Create names for trimmed FASTQ files.

trimming_prefix = os.path.join(fastq_folder, args.sample_name)
trimmed_fastq = trimming_prefix + "_R1_trim.fastq.gz"
trimmed_fastq_R2 = trimming_prefix + "_R2_trim.fastq.gz"
fastqc_folder = os.path.join(outfolder, "fastqc")
fastqc_report = os.path.join(fastqc_folder,
    args.sample_name + "_R1_trim_fastqc.html")
fastqc_report_R2 = os.path.join(fastqc_folder,
    args.sample_name + "_R2_trim_fastqc.html")
if ngstk.check_command(tools.fastqc):
    ngstk.make_dir(fastqc_folder)

pm.info("trimmomatic local_input_files: {}".format(local_input_files))

trim_cmd_chunks = [
    "{trim} {PE} -threads {cores}".format(
        trim=tools.trimmomatic,
        PE="PE" if args.paired_end else "SE",
        cores=pm.cores),
    untrimmed_fastq1,
    untrimmed_fastq2 if args.paired_end else None,
    trimmed_fastq,
    trimming_prefix + "_R1_unpaired.fq" if args.paired_end else None,
    trimmed_fastq_R2 if args.paired_end else None,
    trimming_prefix + "_R2_unpaired.fq" if args.paired_end else None,
    "ILLUMINACLIP:" + adapters + ":2:30:10 MINLEN:36"
]
trim_cmd = build_command(trim_cmd_chunks)

def check_trim():
    pm.info("Evaluating read trimming")

    if args.paired_end and not trimmed_fastq_R2:
        pm.warning("Specified paired-end but no R2 file")

    n_trim = int(ngstk.count_reads(trimmed_fastq, args.paired_end))
    pm.report_result("Trimmed_reads", int(n_trim))
    try:
        rr = int(pm.get_stat("Raw_reads"))
    except:
        pm.warning("Can't calculate trim loss rate without raw read result.")
    else:
        pm.report_result(
            "Trim_loss_rate", round((rr - n_trim) * 100 / rr, 2))

        # Also run a fastqc (if installed/requested)
    if fastqc_folder and os.path.isabs(fastqc_folder):
        ngstk.make_sure_path_exists(fastqc_folder)
    cmd = (tools.fastqc + " --noextract --outdir " +
           fastqc_folder + " " + trimmed_fastq)
    pm.run(cmd, fastqc_report, nofail=False)
    pm.report_object("FastQC report trim r1", fastqc_report)
    if args.paired_end and trimmed_fastq_R2:
        cmd = (tools.fastqc + " --noextract --outdir " +
               fastqc_folder + " " + trimmed_fastq_R2)
        pm.run(cmd, fastqc_report_R2, nofail=False)
        pm.report_object("FastQC report trim r2", fastqc_report_R2)

pm.run(trim_cmd, trimmed_fastq, follow=check_trim)

#clean fastq files
#pm.clean_add(os.path.join(fastq_folder, "*"), conditional=True)
#pm.clean_add(os.path.join(fastq_folder, "*.log"), conditional=True)

# Prepare variables for alignment step
if args.paired_end:
    unmap_fq1 = trimmed_fastq
    unmap_fq2 = trimmed_fastq_R2
else:
    unmap_fq1 = trimmed_fastq
    unmap_fq2 = None

############################################################################
#                          Genome Alignment                                #
############################################################################

pm.timestamp("### Genome Alignment: ")

# Prepare alignment output folder
map_genome_folder = os.path.join(outfolder,
                                 "aligned_" + args.genome_assembly)
ngstk.make_dir(map_genome_folder)

#Prepare QC folder
QC_folder = os.path.join(outfolder, "QC_" + args.genome_assembly)
ngstk.make_dir(QC_folder)

# Primary endpoint file following alignment and deduplication
mapping_genome_bam_star_path = map_genome_folder + "/" + args.sample_name + "."
mapping_genome_bam_star = mapping_genome_bam_star_path + "Aligned.sortedByCoord.out.bam"
mapping_genome_bam_log = mapping_genome_bam_star_path + "Log.final.out"
mapping_genome_bam = mapping_genome_bam_star_path + "bam"
mapping_genome_bam_bw = mapping_genome_bam_star_path + "bw"
mapping_genome_bam_idx = mapping_genome_bam + ".bai"
mapping_genome_bam_dedup = mapping_genome_bam_star_path + "dedup.bam"
mapping_genome_bam_dedup_metrics = mapping_genome_bam_star_path + "dedup.metrics.txt"
mapping_genome_bam_dedup_unique = mapping_genome_bam_star_path + "dedup.unique.bam"
mapping_genome_bam_dedup_unique_idx = mapping_genome_bam_star_path + "dedup.unique.bam.bai"
mapping_genome_bam_dedup_unique_bw = mapping_genome_bam_star_path + "dedup.unique.bw"

#alignment indexes from command line arguments
STAR_RNA_index = args.STAR_RNA_index
STAR_genome_index = args.STAR_genome_index
Bowtie2_index = args.Bowtie2_index

#generate STAR alignment statistics
def check_alignment():
    if mapper =="STAR":
        x = subprocess.check_output("awk -F\| 'NR==6 { print $2 }' " + mapping_genome_bam_log, shell=True)
        ir = int(x.decode().strip())
        pm.report_result("Input_reads", ir)
        x = subprocess.check_output("awk -F\| 'NR==9 { print $2 }' " + mapping_genome_bam_log, shell=True)
        ur = int(x.decode().strip())
        pm.report_result("Unique_reads", ur)
        x = subprocess.check_output("awk -F\| 'NR==24 { print $2 }' " + mapping_genome_bam_log, shell=True)
        mmr = int(x.decode().strip())
        pm.report_result("Multimapped_reads", mmr)
        x = subprocess.check_output("awk -F\| 'NR==26 { print $2 }' " + mapping_genome_bam_log, shell=True)
        mmxr = int(x.decode().strip())
        pm.report_result("Multimapped_too_many_loci_reads", mmxr)
        x1 = subprocess.check_output("awk -F\| 'NR==29 { print $2 }' " + mapping_genome_bam_log, shell=True)
        x2 = subprocess.check_output("awk -F\| 'NR==31 { print $2 }' " + mapping_genome_bam_log, shell=True)
        x3 = subprocess.check_output("awk -F\| 'NR==33 { print $2 }' " + mapping_genome_bam_log, shell=True)
        unr = int(x1.decode().strip()) + int(x2.decode().strip()) + int(x3.decode().strip())
        pm.report_result("Mapped_reads", (ur+mmr+mmxr))
        pm.report_result("Unmapped_reads", unr)
        pm.report_result("Mapping_rate", round((ur+mmr+mmxr)*100/ir,2))
        pm.report_result("Unique_Mapping_rate", round((ur)*100/ir,2))
        pm.report_result("Multi_Mapping_rate", round((mmr+mmxr)*100/ir,2))
        #mitochondrial Reads
        x = subprocess.check_output("samtools view -c " + mapping_genome_bam + " chrM", shell=True)
        chrM = int(x.decode().strip())
        pm.report_result("Mitochondrial_reads", chrM)
        #mitochondrial Reads percentage
        pm.report_result("Mitochondrial_reads_percentage", round(chrM*100/(ur+mmr+mmxr),2))
    else:
        #mapped reads
        x = subprocess.check_output("samtools view -F 4 -c " + mapping_genome_bam, shell=True)
        mapped = int(x.decode().strip())
        pm.report_result("Mapped_reads", mapped)
        #mapping rate
        trimmed = int(pm.get_stat("Trimmed_reads"))
        mr = round(mapped*100/trimmed,2)
        pm.report_result("Mapping_rate", mr)
        #mitochondrial Reads
        x = subprocess.check_output("samtools view -c " + mapping_genome_bam + " chrM", shell=True)
        chrM = int(x.decode().strip())
        pm.report_result("Mitochondrial_reads", chrM)
        #mitochondrial Reads percentage
        pm.report_result("Mitochondrial_reads_percentage", round(chrM*100/trimmed,2))

if args.pipeline_mode == "repeats":
    if args.protocol == "RNA":
        #annotated genomes for RNA-seq
        mapper = "STAR"
        cmd = "STAR" + " --runThreadN " + str(pm.cores)
        cmd += " --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM"
        cmd += " SortedByCoordinate --runMode alignReads --outFilterMultimapNmax 5000"
        cmd += " --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random"
        cmd += " --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --seedSearchStartLmax 30"
        cmd += " --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000"
        cmd += " --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300"
        cmd += " --seedNoneLociPerWindow 1000 --genomeDir " + STAR_RNA_index
        cmd += " --readFilesCommand zcat --readFilesIn " + unmap_fq1
        if args.paired_end:
                    cmd += " " + unmap_fq2 + " "
        cmd += " --outFileNamePrefix " + mapping_genome_bam_star_path
        cmd2 = "mv " + mapping_genome_bam_star + " " + mapping_genome_bam
    else:
        #non-annotated genomes and for others
        mapper = "STAR"
        cmd = "STAR" + " --runThreadN " + str(pm.cores)
        cmd += " --outSAMtype BAM SortedByCoordinate --runMode alignReads --outFilterMultimapNmax 5000"
        cmd += " --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random"
        cmd += " --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --alignIntronMax 1"
        cmd += " --alignMatesGapMax 350 --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000"
        cmd += " --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300"
        cmd += " --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000"
        cmd += " --genomeDir " + STAR_genome_index
        cmd += " --readFilesCommand zcat --readFilesIn " + unmap_fq1
        if args.paired_end:
                    cmd += " " + unmap_fq2 + " "
        cmd += " --outFileNamePrefix " + mapping_genome_bam_star_path
        cmd2 = "mv " + mapping_genome_bam_star + " " + mapping_genome_bam

    cmd3 = tools.samtools + " index " + mapping_genome_bam
    pm.run([cmd, cmd2, cmd3], mapping_genome_bam, follow=check_alignment)
else: #pipeline_mode: genes
    if args.protocol == "RNA": #use STAR for RNA-seq mapping
        #annotated genomes for RNA-seq
        mapper = "STAR"
        cmd = "STAR" + " --runThreadN " + str(pm.cores)
        cmd += " --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM"
        cmd += " SortedByCoordinate --runMode alignReads --genomeDir " + STAR_RNA_index
        cmd += " --readFilesCommand zcat --readFilesIn " + unmap_fq1
        if args.paired_end:
                    cmd += " " + unmap_fq2 + " "
        cmd += " --outFileNamePrefix " + mapping_genome_bam_star_path
        cmd2 = "mv " + mapping_genome_bam_star + " " + mapping_genome_bam
        cmd3 = tools.samtools + " index " + mapping_genome_bam
        pm.run([cmd, cmd2, cmd3], mapping_genome_bam, follow=check_alignment)
    else: #use bowtie2 for gene mapping
        tempdir = tempfile.mkdtemp(dir=map_genome_folder)
        os.chmod(tempdir, 0o771)
        pm.clean_add(tempdir)
        mapper = "Bowtie2"
        cmd = tools.bowtie2 + " -p " + str(pm.cores)
        cmd += " --very-sensitive -X 2000" #mapping options
        cmd += " --rg-id " + args.sample_name
        cmd += " -x " + Bowtie2_index
        if args.paired_end:
            cmd += " -1 " + unmap_fq1 + " -2 " + unmap_fq2
        else:
            cmd += " -U " + unmap_fq1
        cmd += " | " + tools.samtools + " view -bS - -@ 1 "
        cmd += " | " + tools.samtools + " sort - -@ 1"
        cmd += " -T " + tempdir
        cmd += " -o " + mapping_genome_bam
        cmd2 = tools.samtools + " index " + mapping_genome_bam
        pm.run([cmd, cmd2], mapping_genome_bam, follow=check_alignment)

#report BAM file
pm.report_object("BAM_mapped", mapping_genome_bam)

#check insert size distribution
if args.paired_end and (not args.protocol == "RNA"):
    is_prefix = os.path.join(QC_folder, args.sample_name)
    is_file = is_prefix +  ".is.txt"
    is_pdf = is_prefix + ".fraglen.pdf"
    is_png = is_prefix + ".fraglen.png"
    cmd = tools.samtools + " view " + mapping_genome_bam + " | awk '{print $9}' > " + is_file
    pm.run(cmd, is_file)
    #R script for plotting
    fraglen_script = os.path.join(os.path.dirname(__file__), "frag_distribution.R")
    cmd = "Rscript " + fraglen_script + " " + args.sample_name + " " + is_file + " " + QC_folder
    pm.run(cmd, is_pdf)
    pm.report_object("Insert size distribution", is_pdf, anchor_image=is_png)
    #clean file
    pm.clean_add(is_file)

#####################################
#Remove chrM reads for ATAC-seq data
#TODO
#####################################

#Report Mapping statistics
stats = mapping_genome_bam_star_path + "Log.final.out"

def check_duplicates():
    #Report duplication rate (line 8 in metric.txt file)
    if args.paired_end:
        #Read Pair Duplicates (column 7)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $7}' " + mapping_genome_bam_dedup_metrics, shell=True)
        dup = int(x.decode().strip())
        #Read Pair Optical Duplicates (column 8)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $8}' " + mapping_genome_bam_dedup_metrics, shell=True)
        optdup = int(x.decode().strip())
        #Percent Duplication (column 9)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $9}' " + mapping_genome_bam_dedup_metrics, shell=True)
        percdup = float(x.decode().strip())
        #Estimated Library Size (column 10)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $10}' " + mapping_genome_bam_dedup_metrics, shell=True)
        libsize = int(x.decode().strip())
        pm.report_result("Read_pair_duplicates", dup)
        pm.report_result("Read_pair_optical_duplicates", optdup)
        pm.report_result("Percent_duplication", percdup)
        pm.report_result("Estimated_library_size", libsize)
    else:
        #Read Duplicates (column 6)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $6}' " + mapping_genome_bam_dedup_metrics, shell=True)
        dup = int(x.decode().strip())
        #Percent Duplication (column 9)
        x = subprocess.check_output("awk -F'\t' 'NR==8 {print $9}' " + mapping_genome_bam_dedup_metrics, shell=True)
        percdup = float(x.decode().strip())
        pm.report_result("Read_duplicates", dup)
        pm.report_result("Percent_duplication", percdup)

########################################
#clean BAM files if not RNA-seq Analysis
if not args.protocol == "RNA":
    #Remove duplicates
    cmd = tools.picard + " MarkDuplicates --VALIDATION_STRINGENCY LENIENT -I " + mapping_genome_bam
    cmd += " -O " + mapping_genome_bam_dedup + " -M " + mapping_genome_bam_dedup_metrics
    pm.run(cmd, mapping_genome_bam_dedup, follow=check_duplicates)

    #Remove multimapping reads based on MAPQ=255 (STAR) or MAPQ=42 (Bowtie2)
    #for bowtie2 mapping: also remove unmapped reads
    if mapper == "STAR":
        cmd = tools.samtools + " view -b -q 255 " + mapping_genome_bam_dedup + " > " + mapping_genome_bam_dedup_unique
    else: #mapper == Bowtie2
        cmd = tools.samtools + " view -b -F 4 -q 42 " + mapping_genome_bam_dedup + " > " + mapping_genome_bam_dedup_unique
    pm.run(cmd, mapping_genome_bam_dedup_unique)

    #Index deduplicated and unique BAM file
    cmd = tools.samtools + " index " + mapping_genome_bam_dedup_unique
    pm.run(cmd, mapping_genome_bam_dedup_unique_idx)

    #report mapped reads after filtering and filtered mapping rate
    trimmed = int(pm.get_stat("Trimmed_reads"))
    x = subprocess.check_output("samtools view -F 4 -c " + mapping_genome_bam_dedup_unique, shell=True)
    ur = int(x.decode().strip())
    pm.report_result("Mapped_reads_filtered", round(ur,2))
    pm.report_result("Mapping_rate_filtered", round((ur)*100/trimmed,2))

    #Report deduplicated BAM file
    pm.report_object("BAM_dedup_unique", mapping_genome_bam_dedup_unique)

#Generate BigWig files using Deeptools for RNA-seq use mapping_genome_bam file
if args.protocol == "RNA":
    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam
    cmd += " -o " + mapping_genome_bam_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_bw)
    #Report BigWig file
    pm.report_object("BigWig", mapping_genome_bam_bw)
else:
    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam_dedup_unique
    cmd += " -o " + mapping_genome_bam_dedup_unique_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_dedup_unique_bw)
    #Report deduplicated BigWig file
    pm.report_object("BigWig_dedup", mapping_genome_bam_dedup_unique_bw)

#Clean Temporary Files
pm.clean_add(mapping_genome_bam_dedup)

############################################################################
#                    Read filtering for insert size                        #
############################################################################
if (args.protocol == "CT" or args.protocol == "CR") and args.paired_end:
    mapping_genome_bam_dedup_unique_nuc = mapping_genome_bam_star_path + "dedup.unique.nuc.bam"
    mapping_genome_bam_dedup_unique_nuc_idx = mapping_genome_bam_star_path + "dedup.unique.nuc.bam.bai"
    mapping_genome_bam_dedup_unique_nuc_bw = mapping_genome_bam_star_path + "dedup.unique.nuc.bw"
    mapping_genome_bam_dedup_unique_subnuc = mapping_genome_bam_star_path + "dedup.unique.subnuc.bam"
    mapping_genome_bam_dedup_unique_subnuc_idx = mapping_genome_bam_star_path + "dedup.unique.subnuc.bam.bai"
    mapping_genome_bam_dedup_unique_subnuc_bw = mapping_genome_bam_star_path + "dedup.unique.subnuc.bw"

    #nucleosomal reads (insert size >= 120 = insert size^2 > 14400)
    cmd = tools.samtools + " view -h " + mapping_genome_bam_dedup_unique
    cmd += " | awk \'substr($0,1,1)==\"@\" || ($9^2 >= 14400)\' | "
    cmd += tools.samtools + " view -b > " + mapping_genome_bam_dedup_unique_nuc
    pm.run (cmd, mapping_genome_bam_dedup_unique_nuc, shell=True)
    cmd = tools.samtools + " index " + mapping_genome_bam_dedup_unique_nuc
    pm.run (cmd, mapping_genome_bam_dedup_unique_nuc_idx)

    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam_dedup_unique_nuc
    cmd += " -o " + mapping_genome_bam_dedup_unique_nuc_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_dedup_unique_nuc_bw)


    #nucleosomal reads (insert size < 120)
    cmd = tools.samtools + " view -h " + mapping_genome_bam_dedup_unique
    cmd += " | awk \'substr($0,1,1)==\"@\" || ($9^2 < 14400)\' | "
    cmd += tools.samtools + " view -b > " + mapping_genome_bam_dedup_unique_subnuc
    pm.run (cmd, mapping_genome_bam_dedup_unique_subnuc, shell=True)
    cmd = tools.samtools + " index " + mapping_genome_bam_dedup_unique_subnuc
    pm.run (cmd, mapping_genome_bam_dedup_unique_subnuc_idx)

    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam_dedup_unique_subnuc
    cmd += " -o " + mapping_genome_bam_dedup_unique_subnuc_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_dedup_unique_subnuc_bw)


    ############################################################################
    #                          Determine TSS enrichment                        #
    #           taken from PEPATAC Pipeline                                    #
    ############################################################################
if args.protocol == "ATAC":
    if not os.path.exists(args.refgene_tss):
        print("Skipping TSS -- TSS enrichment requires TSS annotation file: {}"
              .format(args.refgene_tss))
    else:
        pm.timestamp("### Calculate TSS enrichment")

        Tss_enrich = os.path.join(QC_folder, args.sample_name +
                                  "_TSS_enrichment.txt")
        cmd = tool_path("pyTssEnrichment.py")
        cmd += " -a " + mapping_genome_bam_dedup_unique + " -b " + args.refgene_tss + " -p ends"
        cmd += " -c " + str(pm.cores)
        cmd += " -z -v -s 6 -o " + Tss_enrich
        print(cmd)
        pm.run(cmd, Tss_enrich, nofail=True)

        if not pm.get_stat('TSS_score') or args.new_start:
            with open(Tss_enrich) as f:
                floats = list(map(float, f))
            try:
                # If the TSS enrichment is 0, don't report
                list_len = 0.05*float(len(floats))
                normTSS = [x / (sum(floats[1:int(list_len)]) /
                           len(floats[1:int(list_len)])) for x in floats]
                max_index = normTSS.index(max(normTSS))

                if (((normTSS[max_index]/normTSS[max_index-1]) > 1.5) and
                    ((normTSS[max_index]/normTSS[max_index+1]) > 1.5)):
                    tmpTSS = list(normTSS)
                    del tmpTSS[max_index]
                    max_index = tmpTSS.index(max(tmpTSS)) + 1

                Tss_score = round(
                    (sum(normTSS[int(max_index-50):int(max_index+50)])) /
                    (len(normTSS[int(max_index-50):int(max_index+50)])), 1)

                pm.report_result("TSS_score", round(Tss_score, 1))
            except ZeroDivisionError:
                pm.report_result("TSS_score", 0)
                pass

        # Call Rscript to plot TSS Enrichment
        Tss_pdf = os.path.join(QC_folder, args.sample_name + "_TSS_enrichment.pdf")
        Tss_png = os.path.join(QC_folder, args.sample_name + "_TSS_enrichment.png")
        cmd = (tools.Rscript + " " + tool_path("plot.TSS.enrichment.R") +
                " -i " + Tss_enrich)
        pm.run(cmd, Tss_pdf, nofail=True)

        pm.report_object("TSS enrichment", Tss_pdf, anchor_image=Tss_png)

#STOP pipeline if genes mode,otherwise continue with repeats coverage
if args.pipeline_mode == "genes":
    pm.stop_pipeline()
    sys.exit()

############################################################################
#                          IAP Coverage                                    #
############################################################################

if args.genome_assembly == "mm10":
    pm.timestamp("### IAP Coverage: ")

    # Prepare output folder
    IAP_coverage_folder = os.path.join(outfolder,
                                 "IAP_coverage")
    ngstk.make_dir(IAP_coverage_folder)

    #processing files
    IAP_plus = IAP_coverage_folder + "/" + args.sample_name + ".IAP.plus.txt"
    IAP_minus = IAP_coverage_folder + "/" + args.sample_name + ".IAP.minus.txt"
    IAP_norm_coverage = IAP_coverage_folder + "/" + args.sample_name + ".IAP.norm.coverage.txt"

    #generate Coverage using bedtools coverage
    #Plus orientation
    cmd = tools.bedtools + " coverage -g " + args.genome_index + " -sorted -d -a " + res.gag_plus
    cmd += " -b " + mapping_genome_bam
    cmd += " > " + IAP_plus
    pm.run(cmd, IAP_plus)

    #Minus Orientation
    cmd = tools.bedtools + " coverage -g " + args.genome_index + " -sorted -d -a " + res.gag_minus
    cmd += " -b " + mapping_genome_bam
    cmd += " > " + IAP_minus
    pm.run(cmd, IAP_minus)

    #summarize and normalize coverage
    norm_factor = float(pm.get_stat("Mapped_reads"))/1000000
    cmd = "tac " + IAP_minus + " | cat " + IAP_plus
    cmd += " | awk -vN=15413 '{s[(NR-1)%N]+=$5}END{for(i=0;i<N;i++){print s[i]/" + str(norm_factor) + "}}'"
    cmd += " > " + IAP_norm_coverage

    pm.run(cmd, IAP_norm_coverage, shell=True)

    pm.clean_add(IAP_plus)
    pm.clean_add(IAP_minus)

############################################################################
#                          Feature Counts                                  #
############################################################################

pm.timestamp("### Feature Counts: ")

# Prepare output folder
feature_counts_folder = os.path.join(outfolder,
                                 "feature_counts")
ngstk.make_dir(feature_counts_folder)

#Target Files
feature_counts_temp = feature_counts_folder + "/" + args.sample_name + ".fc.tmp.txt"
feature_counts_result = feature_counts_folder + "/" + args.sample_name + ".fc.txt"
feature_counts_result_id = feature_counts_folder + "/" + args.sample_name + ".fc.id.txt"

#SAF files
SAF = args.repeats_SAF
SAFid = args.repeats_SAFid #for individual elements

#perform featurecounts on repeat classes
cmd = tools.featureCounts + " -M -F SAF -T 1 -s 0 -a " + SAF
if args.paired_end:
    cmd+= " -p "
cmd += " -o " + feature_counts_temp + " " + mapping_genome_bam

pm.run(cmd, feature_counts_temp)
#reformat output file
cmd = "awk '{print $1,$6,$7}' " + feature_counts_temp + " > " + feature_counts_result
pm.run(cmd, feature_counts_result)
pm.clean_add(feature_counts_temp)

#perform featurecounts on individual repeats (on BAM file with unique mapping)
cmd = tools.featureCounts + " -F SAF -T 1 -s 0 -a " + SAFid
if args.paired_end:
    cmd+= " -p "
cmd += " -o " + feature_counts_result_id + " " + mapping_genome_bam_dedup_unique

pm.run(cmd, feature_counts_result_id)

############################################################################
#                          Pipeline finished                               #
############################################################################

pm.stop_pipeline()
