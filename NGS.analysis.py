#!/usr/bin/env python3
"""
NGS.analysis v3 -- per-sample processing of ChIP-seq, ATAC-seq, CUT&RUN, CUT&Tag and RNA-seq data.

Chromatin (ATAC, CHIP, CR, CT): ONE mode since 3.3.0 (pipeline_mode is ignored).
  bowtie2 (--very-sensitive -X 2000 --dovetail; multimappers at one best locus, chosen at random) -> chrM split
  -> MarkDuplicates -> two fragment sets, both proper pairs, duplicates removed, canonical chromosomes, no chrM:
    filtered (.filt.bam: MAPQ >= FILT_MAPQ = 10, blacklist removed) -> CPM bigwig (also the locus bigwig), peaks +
      FRiP, QC, CR/CT split;
    all (MAPQ >= 0, temporary .all.bam)  -> all-fragment bigwig, repeat family/element counts in two tiers
      (all fragments; locus tier MAPQ >= LOCUS_MAPQ = FILT_MAPQ) with per-element uncertainty flags, IAP coverage (mm10).
  Locus tier and young flag: analysis/NGS.analysis.v3.benchmark/README.md (MAPQ curves, SNP / non-reference TE
  simulation, element calibration); docs/repeats.md.
RNA: genes = STAR + RSEM; repeats = STAR random-1 exactly as in v2 (frozen until the RNA redesign).
All: Trimmomatic with minAdapterLength 1 + keepBothReads (PE), statistics in fragments (docs/stats_units.md).
Chromatin (since 3.1.0): pairs touching the mitochondrial genome are split off straight after alignment, so <s>.bam
and everything derived from it (duplicates, complexity, peaks, FRiP, bigwigs, counts) is nuclear only.
Peaks (since 3.2.0): MACS3 3.1.0 with a numeric effective genome size (genomes/<g>.yaml, MACS3's mm/hs), 01.Angela
settings (broad incl. H3K27ac: -q 1e-4 --broad-cutoff 0.2), ATAC with --call-summits (Peaks_n = regions).
"""

__author__ = ["Gunnar Schotta"]
__email__ = "gunnar.schotta@bmc.med.lmu.de"
__version__ = "3.3.0"

from argparse import ArgumentParser, SUPPRESS
import os
import sys
import tempfile
import subprocess
import logging
from contextlib import contextmanager
import yaml
import pypiper
from pypiper import build_command

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# tools come from the environment of this interpreter (no reliance on an activated env)
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
sys.path.insert(0, SCRIPT_DIR)
import ngs_qc as qc  # noqa: E402

PROTOCOLS = ["CHIP", "CT", "CR", "RNA", "ATAC"]
PIPELINEMODE = ["genes", "repeats"]   # RNA only; chromatin has one mode since 3.3.0
# chromatin MAPQ tiers (bowtie2 gives both mates of a concordant pair the same MAPQ, so a tier is per fragment)
# 10 since 3.3.0 (user decision 2026-11-25, peaks benchmark mapq10/): same peaks/FRiP/replicate concordance as 30 for
# ATAC and H3K27ac, slightly more concordant broad marks; one BAM and one bigwig serve peaks and locus plots.
# Young-copy fragments are 0.1-1.1 % of all fragments of the extra ones. >= 98.5 % correct copy for young L1 at 0.2 %
# strain divergence (2x100); MAPQ 30 adds no precision under divergence.
FILT_MAPQ = 10      # .filt.bam (peaks, FRiP, filtered bigwig, QC)
LOCUS_MAPQ = FILT_MAPQ   # repeat locus tier = the filtered set's cutoff
# young = interspersed repeat < 5 % diverged from its consensus: at MAPQ >= LOCUS_MAPQ these elements receive
# 5-15 % fragments from elsewhere under 129-like divergence, mostly from non-reference copies (element calibration)
YOUNG_DIV_PCT = 5.0
# broad marks; H3K27ac broad since 3.2.0 (01.Angela setting; narrow lost half of the domains, peaks benchmark)
BROAD_TARGETS = {"h3k27me3", "h3k9me3", "h3k36me3", "h3k79me2", "h4k20me3", "h3k9me2", "h2ak119ub", "h3k27ac"}
# MACS3 thresholds (3.2.0 = 01.Angela settings): narrow -q 0.01; broad -q 1e-4 with --broad-cutoff 0.2
MACS_NARROW_ARGS = "-q 0.01"
MACS_BROAD_ARGS = "-q 1e-4 --broad --broad-cutoff 0.2"
NO_PEAK_TARGETS = {"none", "input", "igg", "na", ""}

parser = ArgumentParser(description="NGS analysis version " + __version__)
parser = pypiper.add_pypiper_args(
    parser, groups=["pypiper", "looper", "ngs", "config"],
    required=["input", "genome", "sample-name", "output-parent"])
parser.add_argument("--protocol", dest="protocol", type=str, default=None, choices=PROTOCOLS,
                    help="NGS processing protocol.")
parser.add_argument("--pipeline-mode", dest="pipeline_mode", type=str, default="genes", choices=PIPELINEMODE,
                    help="RNA only: genes (STAR + RSEM) or repeats (v2 STAR random-1). Ignored for chromatin.")
parser.add_argument("--genome-config", dest="genome_config", type=str, default=None,
                    help="Genome resource YAML (default: genomes/<genome>.yaml next to this script).")
parser.add_argument("--STAR_RNA_index", dest="STAR_RNA_index", type=str, default=None, help="STAR RNA index.")
parser.add_argument("--STAR_genome_index", dest="STAR_genome_index", type=str, default=None,
                    help="STAR genome index.")
parser.add_argument("--Bowtie2_index", dest="Bowtie2_index", type=str, default=None, help="Bowtie2 genome index.")
parser.add_argument("--refgene_tss", dest="refgene_tss", type=str, default=None, help="TSS BED file")
parser.add_argument("--genome_index", dest="genome_index", type=str, default=None,
                    help="Genome FASTA index file (.fai)")
parser.add_argument("--repeats_SAF", dest="repeats_SAF", type=str, default=None,
                    help="Repeats SAF file for featureCounts")
parser.add_argument("--repeats_SAFid", dest="repeats_SAFid", type=str, default=None,
                    help="Repeats SAF file for featureCounts, individual repeats")
parser.add_argument("--rsem-index", dest="rsem_index", type=str, default=None,
                    help="RSEM transcriptome index prefix (TPM quantification).")
parser.add_argument("--strandedness", dest="strandedness", type=str, default="reverse",
                    choices=["none", "forward", "reverse"], help="Library strandedness (RSEM).")
parser.add_argument("--peak-mode", dest="peak_mode", type=str, default="auto",
                    choices=["auto", "narrow", "broad", "none"], help="MACS3 peak calling mode.")
parser.add_argument("--target", dest="target", type=str, default=None,
                    help="Antibody target (used by --peak-mode auto).")
parser.add_argument("--se-extend", dest="se_extend", type=int, default=200,
                    help="Fragment length assumed for single-end bigwigs/peaks.")
parser.add_argument("--pipestat-config", dest="pipestat_config", type=str, default=None,
                    help="Looper-generated pipestat config file path.")
parser.add_argument("--legacy-trim", dest="legacy_trim", action="store_true", help=SUPPRESS)  # v2 trimming (validation only)
parser.add_argument("--keep-mito", dest="keep_mito", action="store_true", help=SUPPRESS)  # keep chrM in <s>.bam (v2 equivalence test)

args = parser.parse_args()
args.paired_end = args.single_or_paired.lower() == "paired"
if not args.input:
    parser.print_help()
    raise SystemExit


def _none(x):
    return x in (None, "", "None", "null")


# ----------------------------------------------------------------------------- genome resources
genome_config = args.genome_config or os.path.join(SCRIPT_DIR, "genomes", args.genome_assembly + ".yaml")
GCFG = yaml.safe_load(open(genome_config)) if os.path.exists(genome_config) else {}


def resource(key, argval=None):
    """Explicit command-line value wins over the genome YAML."""
    return argval if not _none(argval) else (None if _none(GCFG.get(key)) else GCFG.get(key))


STAR_RNA_index = resource("star_rna_index", args.STAR_RNA_index)
STAR_genome_index = resource("star_genome_index", args.STAR_genome_index)
Bowtie2_index = resource("bowtie2_index", args.Bowtie2_index)
TSS_BED = resource("tss_bed", args.refgene_tss)
FAI = resource("fai", args.genome_index)
REPEATS_SAF = resource("repeats_saf", args.repeats_SAF)
REPEATS_SAFID = resource("repeats_safid", args.repeats_SAFid)
REPEATS_ANNOT = resource("repeats_annotation")
RSEM_INDEX = resource("rsem_index", args.rsem_index) if args.protocol == "RNA" else None
BLACKLIST = resource("blacklist")
MITO = GCFG.get("mito", "chrM")
CANONICAL = GCFG.get("canonical", r"chr([0-9]+|X|Y)")
MACS_GSIZE = str(GCFG.get("macs_gsize", "hs"))

# ----------------------------------------------------------------------------- output folder, version guard
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))
os.makedirs(outfolder, exist_ok=True)
_marker = os.path.join(outfolder, "NGS.analysis.version")
if os.path.exists(_marker):
    _old = open(_marker).read().strip()
    # major.minor: 3.1 changed <s>.bam (chrM removed), 3.2 the peak calls, 3.3 the chromatin mode (one bowtie2 mode
    # with repeat outputs); pypiper would keep existing targets of an older version
    if _old.split(".")[:2] != __version__.split(".")[:2]:
        sys.exit(f"ERROR: {outfolder} contains results of NGS.analysis {_old}; v{__version__} would silently "
                 f"reuse them (pypiper skips existing targets). Use a new output folder.")
elif any(os.path.exists(os.path.join(outfolder, f)) for f in ("stats.yaml", "NGS.analysis_log.md")):
    sys.exit(f"ERROR: {outfolder} contains results of an older NGS.analysis version (no version marker). "
             f"Use a new output folder.")
with open(_marker, "w") as _f:
    _f.write(__version__ + "\n")

# Per-sample pipestat config: own results file, and ALWAYS the full schema next to this script (the looper-cached
# schema_path can be stale; the old protocol-filtered schema made pipestat reject keys reported in repeats mode).
if args.pipestat_config and os.path.exists(args.pipestat_config):
    with open(args.pipestat_config) as _f:
        _cfg = yaml.safe_load(_f) or {}
    _cfg["results_file_path"] = os.path.join(outfolder, "stats.yaml")
    _cfg["schema_path"] = os.path.join(SCRIPT_DIR, "pipestat_results_schema.yaml")
    _resolved = os.path.join(outfolder, "pipestat_config.yaml")
    with open(_resolved, "w") as _f:
        yaml.dump(_cfg, _f)
    args.pipestat_config = _resolved


@contextmanager
def _skip_duplicate_pypiper_file_handler(log_path):
    """Avoid piper 0.15.1 double-writing logger messages to its markdown log."""
    original_file_handler = logging.FileHandler
    target_log = os.path.abspath(log_path)

    class _NoopPipelineLogHandler(logging.NullHandler):
        baseFilename = target_log

    class _PipelineLogFileHandler(original_file_handler):
        def __new__(cls, filename, *a, **kw):
            if os.path.abspath(filename) == target_log:
                return _NoopPipelineLogHandler()
            return super().__new__(cls)

        def __init__(self, filename, *a, **kw):
            if os.path.abspath(filename) == target_log:
                return
            super().__init__(filename, *a, **kw)

    logging.FileHandler = _PipelineLogFileHandler
    try:
        yield
    finally:
        logging.FileHandler = original_file_handler
        for handler in logging.getLogger().handlers:
            if isinstance(handler, logging.NullHandler) and getattr(handler, "baseFilename", None) == target_log:
                logging.getLogger().removeHandler(handler)
                handler.close()


with _skip_duplicate_pypiper_file_handler(os.path.join(outfolder, "NGS.analysis_log.md")):
    pm = pypiper.PipelineManager(name="NGS.analysis", outfolder=outfolder,
                                 pipestat_record_identifier=args.sample_name,
                                 pipestat_config_file=args.pipestat_config, args=args, version=__version__)
ngstk = pypiper.NGSTk(pm=pm)
tools = pm.config.tools


def tool(name, default):
    v = getattr(tools, name, None) if tools is not None else None
    return v if v else default


SAMTOOLS = tool("samtools", "samtools")
STAR = tool("STAR", "STAR")
RSEM = tool("rsem", "rsem-calculate-expression")
MACS3 = tool("macs3", "macs3")
BEDTOOLS = tool("bedtools", "bedtools")


def tool_path(tool_name):
    return os.path.join(SCRIPT_DIR, tool_name)


def stat(key, cast=float):
    v = pm.get_stat(key)
    return None if v is None else cast(v)


def report(key, value):
    if value is not None:
        pm.report_result(key, value)


# ----------------------------------------------------------------------------- inputs
for f in [args.input[0]] + ([args.input2[0]] if args.input2 else []):
    if not os.path.isfile(f):
        pm.fail_pipeline(IOError("Could not find: {}".format(f)))
    elif os.stat(f).st_size == 0:
        pm.fail_pipeline(IOError("File exists but is empty: {}".format(f)))

report("File_mb", round(ngstk.get_file_size([x for x in [args.input, args.input2] if x is not None]), 2))
report("Read_type", args.single_or_paired)
report("Genome", args.genome_assembly)
report("Protocol", args.protocol)
report("Pipeline_mode", args.pipeline_mode if args.protocol == "RNA" else "chromatin")
report("Pipeline_version", __version__)

raw_folder = os.path.join(outfolder, "raw")
fastq_folder = os.path.join(outfolder, "fastq")

pm.timestamp("### Merge/link and fastq conversion: ")
local_input_files = ngstk.merge_or_link([args.input, args.input2], raw_folder, args.sample_name)
if any(isinstance(i, list) for i in local_input_files):
    local_input_files = [i for e in local_input_files for i in e]
local_input_files = list(dict.fromkeys(local_input_files))
cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(
    local_input_files, args.sample_name, args.paired_end, fastq_folder, zipmode=True)
if any(isinstance(i, list) for i in unaligned_fastq):
    unaligned_fastq = [i for e in unaligned_fastq for i in e]
if any(isinstance(i, dict) for i in local_input_files):
    unaligned_fastq = list(dict.fromkeys(unaligned_fastq))
pm.run(cmd, unaligned_fastq, follow=ngstk.check_fastq(local_input_files, unaligned_fastq, args.paired_end))
pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)
untrimmed_fastq1, untrimmed_fastq2 = (unaligned_fastq[0], unaligned_fastq[1]) if args.paired_end \
    else (unaligned_fastq, None)

fastqc_folder = os.path.join(outfolder, "fastqc")
ngstk.make_sure_path_exists(fastqc_folder)


def fastqc(fq, report_name, html):
    pm.run(tools.fastqc + " --noextract --outdir " + fastqc_folder + " " + fq, html, nofail=False)
    pm.report_object(report_name, html)


fastqc(untrimmed_fastq1, "FastQC report r1", os.path.join(fastqc_folder, args.sample_name + "_R1_fastqc.html"))
if args.paired_end and untrimmed_fastq2:
    fastqc(untrimmed_fastq2, "FastQC report r2", os.path.join(fastqc_folder, args.sample_name + "_R2_fastqc.html"))

############################################################################
#                     Adapter trimming                                     #
############################################################################
pm.timestamp("### Adapter trimming: ")
adapters = tool_path("NexteraPE-PE.fa") if args.protocol in ["ATAC", "CT"] else tool_path("TruSeq3-PE-2.fa")
trimming_prefix = os.path.join(fastq_folder, args.sample_name)
trimmed_fastq = trimming_prefix + "_R1_trim.fastq.gz"
trimmed_fastq_R2 = trimming_prefix + "_R2_trim.fastq.gz"
trim_log = trimming_prefix + ".trimmomatic.log"
if args.legacy_trim:   # v2: palindrome without keepBothReads, minAdapterLength 8 (loses fragments < read length)
    clip, unpaired_ext = ":2:30:10", "_unpaired.fq"
elif args.paired_end:  # v3: clip adapter overhangs down to 1 nt, keep both reads of short fragments
    clip, unpaired_ext = ":2:30:10:1:true", "_unpaired.fq.gz"
else:                  # SE: palindrome mode not available; simple-clip threshold 7 (~12 nt adapter match)
    clip, unpaired_ext = ":2:30:7", "_unpaired.fq.gz"
trim_cmd = build_command([
    "{} {} -threads {}".format(tools.trimmomatic, "PE" if args.paired_end else "SE", pm.cores),
    untrimmed_fastq1,
    untrimmed_fastq2 if args.paired_end else None,
    trimmed_fastq,
    trimming_prefix + "_R1" + unpaired_ext if args.paired_end else None,
    trimmed_fastq_R2 if args.paired_end else None,
    trimming_prefix + "_R2" + unpaired_ext if args.paired_end else None,
    "ILLUMINACLIP:" + adapters + clip + " MINLEN:30"]) + " 2> " + trim_log


def check_trim():
    t = qc.parse_trimmomatic_log(trim_log, args.paired_end)
    report("Trim_input_fragments", t["input"])
    report("Trimmed_fragments", t["both"])
    if args.paired_end:
        report("Trim_R1_only_fragments", t["forward_only"])
        report("Trim_R2_only_fragments", t["reverse_only"])
    report("Trim_dropped_fragments", t["dropped"])
    report("Trimmed_pct", round(100 * t["both"] / t["input"], 2) if t["input"] else 0)
    fastqc(trimmed_fastq, "FastQC report trim r1",
           os.path.join(fastqc_folder, args.sample_name + "_R1_trim_fastqc.html"))
    if args.paired_end:
        fastqc(trimmed_fastq_R2, "FastQC report trim r2",
               os.path.join(fastqc_folder, args.sample_name + "_R2_trim_fastqc.html"))


pm.run(trim_cmd, trimmed_fastq, follow=check_trim, shell=True)
if args.paired_end:
    pm.clean_add(trimming_prefix + "_R1" + unpaired_ext)
    pm.clean_add(trimming_prefix + "_R2" + unpaired_ext)
unmap_fq1, unmap_fq2 = trimmed_fastq, (trimmed_fastq_R2 if args.paired_end else None)

############################################################################
#                          Genome alignment                                #
############################################################################
pm.timestamp("### Genome Alignment: ")
map_genome_folder = os.path.join(outfolder, "aligned_" + args.genome_assembly)
ngstk.make_dir(map_genome_folder)
QC_folder = os.path.join(outfolder, "QC_" + args.genome_assembly)
ngstk.make_dir(QC_folder)

prefix = map_genome_folder + "/" + args.sample_name + "."
mapping_genome_bam_star = prefix + "Aligned.out.bam"
mapping_genome_bam_log = prefix + "Log.final.out"
mapping_genome_bam = prefix + "bam"
mapping_genome_bam_bw = prefix + "bw"
mapping_genome_bam_dedup = prefix + "dedup.bam"
mapping_genome_bam_dedup_metrics = prefix + "dedup.metrics.txt"
mapping_genome_bam_dedup_unique = prefix + "dedup.unique.bam"      # RNA repeats mode (v2 name and content)
mapping_genome_bam_dedup_unique_idx = prefix + "dedup.unique.bam.bai"
mapping_genome_bam_dedup_unique_bw = prefix + "dedup.unique.bw"
filt_bam = prefix + "filt.bam"                                     # chromatin: MAPQ >= FILT_MAPQ
filt_bw = prefix + "filt.bw"
all_bam = prefix + "all.bam"                                       # chromatin: all MAPQ (temporary)
all_bw = prefix + "all.bw"
mito_bam = prefix + MITO + ".bam"                                  # split-off mitochondrial fragments (temporary)

CHROMATIN = args.protocol != "RNA"
RNA_REPEATS = args.protocol == "RNA" and args.pipeline_mode == "repeats"
mapper = "Bowtie2" if CHROMATIN else "STAR"
UNIQUE_MAPQ = FILT_MAPQ if CHROMATIN else 255
# Chromatin: remove fragments touching chrM straight after alignment (aligned against the full genome, so MAPQ
# arbitrates nuclear NUMT vs chrM), before <s>.bam is written. Duplication, complexity, peaks, FRiP, bigwigs and
# counts are then nuclear only. RNA keeps chrM (mitochondrial transcripts are signal).
SPLIT_MITO = (args.protocol != "RNA" and not args.keep_mito and bool(FAI) and os.path.exists(FAI)
              and qc.chrom_present(FAI, MITO))
if SPLIT_MITO:
    pm.clean_add(mito_bam)


def check_alignment():
    if mapper == "STAR":
        s = qc.parse_star_log(mapping_genome_bam_log)
        aligned = s["unique"] + s["multi"]          # "too many loci" reads are not in the BAM -> not aligned
        report("Aligned_fragments", aligned)
        report("Unique_fragments", s["unique"])
        report("Multimapped_fragments", s["multi"])
        report("Too_many_loci_fragments", s["too_many_loci"])
        report("Unmapped_fragments", s["input"] - aligned)
        report("Alignment_pct", round(100 * aligned / s["input"], 2) if s["input"] else 0)
        report("Unique_pct", round(100 * s["unique"] / s["input"], 2) if s["input"] else 0)
        report("Multimapped_pct", round(100 * s["multi"] / s["input"], 2) if s["input"] else 0)
    else:
        # bowtie2: counted from the BAMs (with SPLIT_MITO: nuclear <s>.bam + split-off mito BAM)
        bams = [mapping_genome_bam] + ([mito_bam] if SPLIT_MITO else [])
        trimmed = stat("Trimmed_fragments", int)
        aligned = sum(qc.count_fragments(SAMTOOLS, b, args.paired_end, flags_excl=2308) for b in bams)
        uniq = sum(qc.count_fragments(SAMTOOLS, b, args.paired_end, flags_excl=2308, mapq=FILT_MAPQ) for b in bams)
        report("Aligned_fragments", aligned)
        report("Unique_fragments", uniq)
        report("Alignment_pct", round(100 * aligned / trimmed, 2) if trimmed else 0)
        report("Unique_pct", round(100 * uniq / trimmed, 2) if trimmed else 0)
        if args.paired_end:
            report("Aligned_proper_pairs", sum(qc.count_fragments(SAMTOOLS, b, True, flags_req=2, flags_excl=2308)
                                               for b in bams))
    if SPLIT_MITO:
        mito = qc.count_fragments(SAMTOOLS, mito_bam, args.paired_end, flags_excl=2308)
        report("Mito_filter", "removed")
        report("Nuclear_fragments", qc.count_fragments(SAMTOOLS, mapping_genome_bam, args.paired_end, flags_excl=2308))
        report("Mito_fragments", mito)
        report("Mito_pct", round(100 * mito / aligned, 2) if aligned else 0)
    elif FAI and qc.chrom_present(FAI, MITO):
        mito = qc.count_fragments(SAMTOOLS, mapping_genome_bam, args.paired_end, flags_excl=2308, region=MITO)
        report("Mito_filter", "kept")
        report("Mito_fragments", mito)
        report("Mito_pct", round(100 * mito / aligned, 2) if aligned else 0)


if RNA_REPEATS:
    # ---- v2 STAR RNA command, unchanged (Teissandier 2019 strategy; frozen until the RNA redesign) ----
    cmd = STAR + " --runThreadN " + str(pm.cores)
    cmd += " --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM"
    cmd += " Unsorted --runMode alignReads --outFilterMultimapNmax 5000"
    cmd += " --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random"
    cmd += " --winAnchorMultimapNmax 5000 --alignEndsType EndToEnd --seedSearchStartLmax 30"
    cmd += " --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000"
    cmd += " --alignTranscriptsPerWindowNmax 300 --seedPerReadNmax 3000 --seedPerWindowNmax 300"
    cmd += " --outSAMattrRGline ID:" + args.sample_name + " SM:" + args.sample_name
    cmd += " --seedNoneLociPerWindow 1000 --genomeDir " + STAR_RNA_index
    cmd += " --readFilesCommand zcat --readFilesIn " + unmap_fq1
    if args.paired_end:
        cmd += " " + unmap_fq2 + " "
    cmd += " --outFileNamePrefix " + prefix
    cmd2 = SAMTOOLS + " sort " + mapping_genome_bam_star + " -o " + mapping_genome_bam
    cmd3 = SAMTOOLS + " index " + mapping_genome_bam
    pm.run([cmd, cmd2, cmd3], mapping_genome_bam, follow=check_alignment)
elif args.protocol == "RNA":
    cmd = STAR + " --runThreadN " + str(pm.cores)
    cmd += " --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM"
    cmd += " Unsorted --runMode alignReads --genomeDir " + STAR_RNA_index
    cmd += " --readFilesCommand zcat --readFilesIn " + unmap_fq1
    if args.paired_end:
        cmd += " " + unmap_fq2 + " "
    cmd += " --outFileNamePrefix " + prefix
    cmd2 = SAMTOOLS + " sort " + mapping_genome_bam_star + " -o " + mapping_genome_bam
    cmd3 = SAMTOOLS + " index " + mapping_genome_bam
    pm.run([cmd, cmd2, cmd3], mapping_genome_bam, follow=check_alignment)
else:
    tempdir = tempfile.mkdtemp(dir=map_genome_folder)
    os.chmod(tempdir, 0o771)
    pm.clean_add(tempdir)
    cmd = tools.bowtie2 + " -p " + str(pm.cores)
    cmd += " --very-sensitive -X 2000 --dovetail"  # --dovetail: fully clipped short pairs stay concordant
    cmd += " --rg-id " + args.sample_name + " --rg SM:" + args.sample_name
    cmd += " -x " + Bowtie2_index
    cmd += (" -1 " + unmap_fq1 + " -2 " + unmap_fq2) if args.paired_end else (" -U " + unmap_fq1)
    cmd += " 2> " + prefix + "bowtie2.log"
    if SPLIT_MITO:
        cmd += " | " + qc.mito_filter_cmd(SAMTOOLS, MITO, mito_bam)
    cmd += " | " + SAMTOOLS + " sort -@ " + str(max(1, pm.cores // 4)) + " -m 2G -T " + tempdir + "/sort -o "
    cmd += mapping_genome_bam + " -"
    cmd2 = SAMTOOLS + " index " + mapping_genome_bam
    pm.run([cmd, cmd2], mapping_genome_bam, follow=check_alignment, shell=True)

pm.report_object("BAM_mapped", mapping_genome_bam)

############################################################################
#                     RSEM quantification (RNA-seq)                        #
############################################################################
if args.protocol == "RNA" and RSEM_INDEX:
    pm.timestamp("### RSEM quantification: ")
    transcriptome_bam = prefix + "Aligned.toTranscriptome.out.bam"
    rsem_output_prefix = os.path.join(map_genome_folder, args.sample_name + ".rsem")
    rsem_genes_results = rsem_output_prefix + ".genes.results"
    rsem_isoforms_results = rsem_output_prefix + ".isoforms.results"
    cmd_rsem = RSEM + " --alignments"
    if args.paired_end:
        cmd_rsem += " --paired-end"
    cmd_rsem += " --strandedness " + args.strandedness + " --no-bam-output -p " + str(pm.cores)
    cmd_rsem += " " + transcriptome_bam + " " + RSEM_INDEX + " " + rsem_output_prefix
    pm.run(cmd_rsem, rsem_genes_results)
    pm.report_object("RSEM_genes_results", rsem_genes_results)
    pm.report_object("RSEM_isoforms_results", rsem_isoforms_results)

############################################################################
#                Duplicates, filtering, bigwigs                            #
############################################################################


def check_duplicates():
    d = qc.parse_picard_dup_metrics(mapping_genome_bam_dedup_metrics)
    if args.paired_end:
        report("Duplicate_fragments", d["READ_PAIR_DUPLICATES"])
        report("Optical_duplicate_fragments", d["READ_PAIR_OPTICAL_DUPLICATES"])
    else:
        report("Duplicate_fragments", d["UNPAIRED_READ_DUPLICATES"])
    report("Duplication_pct", round(d["PERCENT_DUPLICATION"], 2))


analysis_bam = None  # the BAM used for QC/peaks downstream
if CHROMATIN or RNA_REPEATS:
    cmd = tools.picard + " MarkDuplicates --VALIDATION_STRINGENCY LENIENT -I " + mapping_genome_bam
    cmd += " -O " + mapping_genome_bam_dedup + " -M " + mapping_genome_bam_dedup_metrics
    pm.run(cmd, mapping_genome_bam_dedup, follow=check_duplicates)

if RNA_REPEATS:
    # ---- v2 filter, unchanged: unique (MAPQ 255) records; duplicates are MARKED, not removed ----
    cmd = SAMTOOLS + " view -b -q 255 " + mapping_genome_bam_dedup + " > " + mapping_genome_bam_dedup_unique
    pm.run(cmd, mapping_genome_bam_dedup_unique)
    pm.run(SAMTOOLS + " index " + mapping_genome_bam_dedup_unique, mapping_genome_bam_dedup_unique_idx)
    analysis_bam = mapping_genome_bam_dedup_unique
    report("Filtered_fragments", qc.count_fragments(SAMTOOLS, analysis_bam, args.paired_end, flags_excl=2308))
    report("Filtered_nondup_fragments",
           qc.count_fragments(SAMTOOLS, analysis_bam, args.paired_end, flags_excl=3332))
    pm.report_object("BAM_dedup_unique", mapping_genome_bam_dedup_unique)
elif CHROMATIN:
    # ---- ENCODE-style filtered BAM (MAPQ >= FILT_MAPQ, blacklist removed) and the all-fragment BAM (MAPQ >= 0,
    # temporary, for repeat counts and the all-fragment bigwig). Both: proper pairs, no duplicates, canonical, no chrM.
    canon = qc.canonical_bed(FAI, CANONICAL, MITO, os.path.join(QC_folder, "canonical_chroms.bed"))
    n = str(pm.cores)
    ftmp = tempfile.mkdtemp(dir=map_genome_folder)
    pm.clean_add(ftmp)
    bl = (" | " + BEDTOOLS + " intersect -v -abam stdin -b " + BLACKLIST) if BLACKLIST else ""
    if args.paired_end:
        cmd = (SAMTOOLS + " view -u -f 2 -F 1804 -q " + str(FILT_MAPQ) + " -L " + canon + " " +
               mapping_genome_bam_dedup + bl +
               " | " + SAMTOOLS + " sort -n -@ " + n + " -m 1G -T " + ftmp + "/n -O BAM -" +
               " | " + SAMTOOLS + " fixmate -r - -" +                        # drop mates of removed reads
               " | " + SAMTOOLS + " view -u -f 2 -F 1804 -" +
               " | " + SAMTOOLS + " sort -@ " + n + " -m 1G -T " + ftmp + "/c -o " + filt_bam + " -")
        # flags are set on both mates (proper pair, duplicate), so no orphan mates: the coordinate order is kept
        cmd_all = SAMTOOLS + " view -b -f 2 -F 1804 -L " + canon + " -o " + all_bam + " " + mapping_genome_bam_dedup
    else:
        cmd = (SAMTOOLS + " view -u -F 1804 -q " + str(FILT_MAPQ) + " -L " + canon + " " + mapping_genome_bam_dedup +
               bl + " | " + SAMTOOLS + " sort -@ " + n + " -m 1G -T " + ftmp + "/c -o " + filt_bam + " -")
        cmd_all = SAMTOOLS + " view -b -F 1804 -L " + canon + " -o " + all_bam + " " + mapping_genome_bam_dedup
    pm.run([cmd, SAMTOOLS + " index " + filt_bam], filt_bam, shell=True)
    analysis_bam = filt_bam
    nf = qc.count_fragments(SAMTOOLS, filt_bam, args.paired_end)
    report("Filtered_fragments", nf)
    tr = stat("Trimmed_fragments", int)
    report("Filtered_pct", round(100 * nf / tr, 2) if tr else 0)
    pm.report_object("BAM_filtered", filt_bam)
    pm.run([cmd_all, SAMTOOLS + " index " + all_bam], all_bam + ".bai")
    pm.clean_add(all_bam)
    pm.clean_add(all_bam + ".bai")
    report("All_fragments", qc.count_fragments(SAMTOOLS, all_bam, args.paired_end))

if args.protocol == "RNA" and not RNA_REPEATS:
    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam
    cmd += " -o " + mapping_genome_bam_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_bw)
    pm.report_object("BigWig", mapping_genome_bam_bw)
elif RNA_REPEATS:
    cmd = tools.bamcoverage + " --bam " + mapping_genome_bam_dedup_unique
    cmd += " -o " + mapping_genome_bam_dedup_unique_bw + " --binSize 10 --normalizeUsing RPKM"
    pm.run(cmd, mapping_genome_bam_dedup_unique_bw)
    pm.report_object("BigWig_dedup", mapping_genome_bam_dedup_unique_bw)
else:
    # CPM, fragment-extended, each tier normalised to its own fragment count
    ext = " --extendReads" if args.paired_end else " --extendReads " + str(args.se_extend)
    bwc = " --binSize 10 --normalizeUsing CPM" + ext + " -p " + str(pm.cores)
    pm.run(tools.bamcoverage + " --bam " + filt_bam + " -o " + filt_bw + bwc, filt_bw)
    pm.report_object("BigWig_filtered", filt_bw)
    pm.run(tools.bamcoverage + " --bam " + all_bam + " -o " + all_bw + bwc, all_bw)
    pm.report_object("BigWig_all", all_bw)

pm.clean_add(mapping_genome_bam_dedup)

############################################################################
#                CUT&RUN / CUT&Tag: (sub)nucleosomal split                 #
############################################################################
if (args.protocol == "CT" or args.protocol == "CR") and args.paired_end:
    src, tag, norm, extra = filt_bam, "filt", "CPM", " --extendReads -p " + str(pm.cores)
    cond_sub = "($9^2 < 14400 && $9 != 0)"          # TLEN 0 records go to neither part
    for part, cond in [("nuc", "($9^2 >= 14400)"), ("subnuc", cond_sub)]:
        pbam, pbw = prefix + tag + "." + part + ".bam", prefix + tag + "." + part + ".bw"
        cmd = SAMTOOLS + " view -h " + src + " | awk \'substr($0,1,1)==\"@\" || " + cond + "\' | "
        cmd += SAMTOOLS + " view -b > " + pbam
        pm.run(cmd, pbam, shell=True)
        pm.run(SAMTOOLS + " index " + pbam, pbam + ".bai")
        pm.run(tools.bamcoverage + " --bam " + pbam + " -o " + pbw + " --binSize 10 --normalizeUsing " + norm + extra,
               pbw)
        pm.report_object("BigWig_filtered_" + part, pbw)

############################################################################
#                QC: insert size, library complexity, TSS                  #
############################################################################
if args.protocol != "RNA" and analysis_bam:
    if args.paired_end and stat("Insert_size_median") is None:
        is_prefix = os.path.join(QC_folder, args.sample_name)
        for k, v in qc.insert_size_qc(SAMTOOLS, analysis_bam, is_prefix).items():
            report(k, v)
        if os.path.exists(is_prefix + ".fraglen.pdf"):
            pm.report_object("Insert size distribution", is_prefix + ".fraglen.pdf",
                             anchor_image=is_prefix + ".fraglen.png")
    if stat("NRF") is None:
        ctmp = tempfile.mkdtemp(dir=QC_folder)
        for k, v in qc.library_complexity(mapping_genome_bam, args.paired_end, UNIQUE_MAPQ, ctmp,
                                          threads=pm.cores).items():
            report(k, v)
        subprocess.call(["rm", "-rf", ctmp])

if args.protocol == "ATAC" and analysis_bam:
    if not TSS_BED or not os.path.exists(TSS_BED):
        pm.info("Skipping TSS enrichment: no TSS annotation for genome " + args.genome_assembly)
    else:
        pm.timestamp("### Calculate TSS enrichment")
        tss_prefix = os.path.join(QC_folder, args.sample_name)
        tss_profile = tss_prefix + "_TSS_enrichment.txt"
        cmd = sys.executable + " " + tool_path("pyTssEnrichment.py") + " -a " + analysis_bam + " -b " + TSS_BED
        cmd += " -p ends -c " + str(pm.cores) + " -z -v -s 6 -q " + str(min(UNIQUE_MAPQ, 30)) + " -o " + tss_profile
        pm.run(cmd, tss_profile, nofail=True)
        if os.path.exists(tss_profile):
            score, norm = qc.tss_score(tss_profile)
            report("TSS_score", score)
            if norm is not None:
                qc.plot_tss(norm, score, args.sample_name, tss_prefix)
                pm.report_object("TSS enrichment", tss_prefix + "_TSS_enrichment.pdf",
                                 anchor_image=tss_prefix + "_TSS_enrichment.png")

############################################################################
#                Peaks + FRiP (no control; controls: NGS.peaks.py)         #
############################################################################
peak_mode = args.peak_mode
if peak_mode == "auto":
    tgt = (args.target or "").strip().lower()
    peak_mode = ("none" if tgt in NO_PEAK_TARGETS and args.protocol != "ATAC"
                 else "broad" if tgt in BROAD_TARGETS else "narrow")
if args.protocol != "RNA" and analysis_bam:
    report("Peak_mode", peak_mode)
if args.protocol != "RNA" and analysis_bam and peak_mode != "none":
    pm.timestamp("### Peak calling (MACS3, " + peak_mode + ")")
    peak_folder = os.path.join(outfolder, "peaks_" + args.genome_assembly)
    ngstk.make_dir(peak_folder)
    ptype = "broadPeak" if peak_mode == "broad" else "narrowPeak"
    peaks = os.path.join(peak_folder, args.sample_name + "_peaks." + ptype)
    peaks_nobl = os.path.join(peak_folder, args.sample_name + "_peaks_noBL." + ptype)
    keepdup = "all"                                     # duplicates already removed from .filt.bam
    cmd = MACS3 + " callpeak -t " + analysis_bam + " -n " + args.sample_name + " --outdir " + peak_folder
    cmd += " -g " + MACS_GSIZE + " --keep-dup " + keepdup
    cmd += " -f BAMPE" if args.paired_end else " -f BAM --nomodel --extsize " + str(args.se_extend)
    cmd += " " + (MACS_BROAD_ARGS if peak_mode == "broad" else MACS_NARROW_ARGS)
    # ATAC: one record per summit, so multi-summit regions can be split into individual peaks downstream
    # (summit-based consensus, e.g. make_common_peakset); region coordinates are repeated per summit
    call_summits = args.protocol == "ATAC" and peak_mode == "narrow"
    if call_summits:
        cmd += " --call-summits"
    cmd += " 2> " + os.path.join(peak_folder, args.sample_name + ".macs3.log")
    pm.run(cmd, peaks, shell=True, nofail=True)
    if os.path.exists(peaks):
        cmd = (BEDTOOLS + " intersect -v -a " + peaks + " -b " + BLACKLIST + " > " + peaks_nobl) if BLACKLIST \
            else ("cp " + peaks + " " + peaks_nobl)
        pm.run(cmd, peaks_nobl, shell=True)
        n_peaks = len(qc.peak_regions(peaks_nobl))          # regions, not summit records
        report("Peaks_n", n_peaks)
        if call_summits:
            report("Peak_summits_n", qc.count_lines(peaks_nobl))
    if os.path.exists(peaks) and n_peaks == 0:
        pm.info("No peaks called; FRiP not computed")
    elif os.path.exists(peaks):
        saf = os.path.join(peak_folder, args.sample_name + "_peaks_noBL.saf")
        qc.peaks_to_saf(peaks_nobl, saf)
        frip_out = os.path.join(peak_folder, args.sample_name + ".frip.txt")
        cmd = tools.featureCounts + " -F SAF -a " + saf + " -o " + frip_out + " -T " + str(pm.cores) + " --ignoreDup"
        cmd += (" -p --countReadPairs " if args.paired_end else " ") + analysis_bam
        pm.run(cmd, frip_out + ".summary", nofail=True)
        if os.path.exists(frip_out + ".summary"):
            s = qc.parse_featurecounts_summary(frip_out + ".summary")
            usable = sum(s.get(k, 0) for k in ["Assigned", "Unassigned_NoFeatures", "Unassigned_Ambiguity",
                                               "Unassigned_Overlapping_Length"])
            report("FRiP_pct", round(100 * s.get("Assigned", 0) / usable, 2) if usable else 0)
        pm.report_object("Peaks", peaks_nobl)

############################################################################
#        Chromatin repeats: family / element counts, flags, IAP coverage   #
############################################################################
if CHROMATIN and REPEATS_SAF and REPEATS_SAFID:
    pm.timestamp("### Repeat counts (all fragments / locus tier MAPQ >= " + str(LOCUS_MAPQ) + ")")
    rep_folder = os.path.join(outfolder, "repeats_" + args.genome_assembly)
    ngstk.make_dir(rep_folder)
    rp = os.path.join(rep_folder, args.sample_name + ".")
    # fragments (PE: --countReadPairs) of the all-fragment BAM; -O: a fragment overlapping several annotation entries
    # (nested / overlapping) counts for each of them; locus tier: -Q (concordant mates share their MAPQ)
    fc = (tools.featureCounts + " -F SAF -O -T " + str(pm.cores) +
          (" -p --countReadPairs -B -C" if args.paired_end else ""))
    out_fc = {}
    for level, saf in (("family", REPEATS_SAF), ("element", REPEATS_SAFID)):
        for tier, q in (("all", ""), ("locus", " -Q " + str(LOCUS_MAPQ))):
            f = rp + level + "." + tier + ".fc.txt"
            pm.run(fc + q + " -a " + saf + " -o " + f + " " + all_bam, f + ".summary")
            out_fc[(level, tier)] = f
    elements_tsv = rp + "repeats.elements.tsv.gz"
    family_tsv = rp + "repeats.family.tsv"
    annot = REPEATS_ANNOT if REPEATS_ANNOT and os.path.exists(REPEATS_ANNOT) else None
    if not annot:
        pm.info("No repeats_annotation for " + args.genome_assembly + ": element table without divergence/flags")
    if not os.path.exists(elements_tsv):
        pm.timestamp("### Repeat tables (+ annotation, flags)")
        rsum = qc.repeat_tables(out_fc[("element", "all")], out_fc[("element", "locus")], out_fc[("family", "all")],
                                out_fc[("family", "locus")], annot, elements_tsv + ".tmp.gz", family_tsv,
                                YOUNG_DIV_PCT)
        os.replace(elements_tsv + ".tmp.gz", elements_tsv)
        if annot:
            tot, young = rsum["element_locus_assignments"], rsum["element_young_locus_assignments"]
            report("Young_locus_pct", round(100 * young / tot, 2) if tot else 0)
    fa = qc.parse_featurecounts_summary(out_fc[("family", "all")] + ".summary")
    fl = qc.parse_featurecounts_summary(out_fc[("family", "locus")] + ".summary")
    nall = stat("All_fragments", int)
    report("Repeat_fragments", fa.get("Assigned", 0))
    report("Repeat_pct", round(100 * fa.get("Assigned", 0) / nall, 2) if nall else 0)
    report("Repeat_locus_fragments", fl.get("Assigned", 0))
    report("Repeat_locus_pct", round(100 * fl.get("Assigned", 0) / fa["Assigned"], 2) if fa.get("Assigned") else 0)
    pm.report_object("Repeat_family_table", family_tsv)
    pm.report_object("Repeat_element_table", elements_tsv)
    for f in out_fc.values():
        pm.clean_add(f)          # merged into the two tables (the .summary files are kept)

    # IAP gag coverage (mm10 resource): per-base read depth of the all-fragment BAM (every mate counts, as in v2),
    # per million all-fragment-BAM fragments (v2: per million STAR unique + multi + too-many-loci pairs, incl. chrM)
    if GCFG.get("iap_plus_bed") and nall:
        pm.timestamp("### IAP Coverage: ")
        iap_folder = os.path.join(outfolder, "IAP_coverage")
        ngstk.make_dir(iap_folder)
        IAP_plus = os.path.join(iap_folder, args.sample_name + ".IAP.plus.txt")
        IAP_minus = os.path.join(iap_folder, args.sample_name + ".IAP.minus.txt")
        IAP_norm_coverage = os.path.join(iap_folder, args.sample_name + ".IAP.norm.coverage.txt")
        for bed, out in ((GCFG["iap_plus_bed"], IAP_plus), (GCFG["iap_minus_bed"], IAP_minus)):
            pm.run(BEDTOOLS + " coverage -g " + FAI + " -sorted -d -a " + tool_path(bed) + " -b " + all_bam +
                   " > " + out, out, shell=True)
        cmd = "tac " + IAP_minus + " | cat " + IAP_plus
        cmd += " | awk -vN=15413 '{s[(NR-1)%N]+=$5}END{for(i=0;i<N;i++){print s[i]/" + str(nall / 1e6) + "}}'"
        cmd += " > " + IAP_norm_coverage
        pm.run(cmd, IAP_norm_coverage, shell=True)
        pm.report_object("IAP_coverage", IAP_norm_coverage)
        pm.clean_add(IAP_plus)
        pm.clean_add(IAP_minus)

if not RNA_REPEATS:
    pm.stop_pipeline()
    sys.exit()

############################################################################
#             RNA repeats mode: IAP Coverage (v2, unchanged)               #
############################################################################
if args.genome_assembly == "mm10" and GCFG.get("iap_plus_bed"):
    pm.timestamp("### IAP Coverage: ")
    IAP_coverage_folder = os.path.join(outfolder, "IAP_coverage")
    ngstk.make_dir(IAP_coverage_folder)
    IAP_plus = IAP_coverage_folder + "/" + args.sample_name + ".IAP.plus.txt"
    IAP_minus = IAP_coverage_folder + "/" + args.sample_name + ".IAP.minus.txt"
    IAP_norm_coverage = IAP_coverage_folder + "/" + args.sample_name + ".IAP.norm.coverage.txt"
    cmd = tools.bedtools + " coverage -g " + FAI + " -sorted -d -a " + tool_path(GCFG["iap_plus_bed"])
    cmd += " -b " + mapping_genome_bam + " > " + IAP_plus
    pm.run(cmd, IAP_plus)
    cmd = tools.bedtools + " coverage -g " + FAI + " -sorted -d -a " + tool_path(GCFG["iap_minus_bed"])
    cmd += " -b " + mapping_genome_bam + " > " + IAP_minus
    pm.run(cmd, IAP_minus)
    # v2 normalisation (frozen): STAR unique + multi + too-many-loci (pairs for PE), per million
    _s = qc.parse_star_log(mapping_genome_bam_log)
    norm_factor = float(_s["unique"] + _s["multi"] + _s["too_many_loci"]) / 1000000
    cmd = "tac " + IAP_minus + " | cat " + IAP_plus
    cmd += " | awk -vN=15413 '{s[(NR-1)%N]+=$5}END{for(i=0;i<N;i++){print s[i]/" + str(norm_factor) + "}}'"
    cmd += " > " + IAP_norm_coverage
    pm.run(cmd, IAP_norm_coverage, shell=True)
    pm.clean_add(IAP_plus)
    pm.clean_add(IAP_minus)

############################################################################
#             RNA repeats mode: Feature Counts (v2, unchanged)             #
############################################################################
pm.timestamp("### Feature Counts: ")
feature_counts_folder = os.path.join(outfolder, "feature_counts")
ngstk.make_dir(feature_counts_folder)
feature_counts_temp = feature_counts_folder + "/" + args.sample_name + ".fc.tmp.txt"
feature_counts_result = feature_counts_folder + "/" + args.sample_name + ".fc.txt"
feature_counts_result_id = feature_counts_folder + "/" + args.sample_name + ".fc.id.txt"
cmd = tools.featureCounts + " -M -F SAF -T 1 -s 0 -a " + REPEATS_SAF
if args.paired_end:
    cmd += " -p "
cmd += " -o " + feature_counts_temp + " " + mapping_genome_bam
pm.run(cmd, feature_counts_temp)
cmd = "awk '{print $1,$6,$7}' " + feature_counts_temp + " > " + feature_counts_result
pm.run(cmd, feature_counts_result)
pm.clean_add(feature_counts_temp)
cmd = tools.featureCounts + " -F SAF -T 1 -s 0 -a " + REPEATS_SAFID
if args.paired_end:
    cmd += " -p "
cmd += " -o " + feature_counts_result_id + " " + mapping_genome_bam_dedup_unique
pm.run(cmd, feature_counts_result_id)

pm.stop_pipeline()
