#! /usr/bin/env Rscript
#
# NGS.summarizer.R
#
# Interface to produce project level summary files and reports
# for NGS output when called using `looper`
# usage: Rscript /path/to/Rscript/NGS.summarizer.R
#        /path/to/project_config.yaml
# Depends: R (>= 3.5.1)
# Imports: PEPATACr, argparser
###############################################################################
##### LOAD ARGUMENTPARSER #####
loadLibrary <- tryCatch (
    {
        suppressWarnings(suppressPackageStartupMessages(library(argparser)))
    },
    error=function(e) {
        message("Error: Install the \"argparser\"",
                " library before proceeding.")
        return(NULL)
    },
    warning=function(e) {
        message(e)
        return(TRUE)
    }
)
if (length(loadLibrary)!=0) {
    suppressWarnings(library(argparser))
} else {
    quit()
}
library(pepr)
library(data.table)
library(reshape2)
library(ggplot2)
library(yaml) 
library(SummarizedExperiment)

# Create a parser
p <- arg_parser("Produce Summary Reports, Files, and Plots")

# Add command line arguments
p <- add_argument(p, arg="config", short="-c",
                  help="project_config.yaml")
p <- add_argument(p, arg="output", short="-o",
                  help="Project parent output directory path")
p <- add_argument(p, arg="results", short="-r",
                  help="Project results output subdirectory path")
p <- add_argument(p, arg="--new-start", short="-N", flag=TRUE,
                  help=paste0("New start mode. This flag will tell the ",
                       "summarizer to start over, and run every command, even ",
                       "if its target output already exists."))

# Parse the command line arguments
argv <- parse_args(p)

#' Create and return sample statistics summary data table
#'
#' @param samples A PEP project character vector of sample names
#' @param results_subdir A PEP project results subdirectory path
#' @export

createStatsSummary <- function(samples, results_subdir) {
    write(paste0("Creating stats summary from YAML..."), stdout())
    stats <- NULL
    missing_files <- 0

    for (sample in samples) {
        sample_dir  <- file.path(results_subdir, sample)
        yaml_file   <- file.path(sample_dir, "stats.yaml")
        if (!file.exists(yaml_file)) {
            missing_files <- missing_files + 1
            next
        }
        y <- yaml::read_yaml(yaml_file)
        # Assume a single top-level pipeline key (e.g., "NGS.analysis")
        pipe_name <- names(y)[1]
        node <- tryCatch(y[[pipe_name]][["sample"]][[sample]], error = function(e) NULL)
        if (is.null(node)) {
            warning(sprintf("No stats for sample %s in %s", sample, yaml_file))
            next
        }
        # Flatten one sample's list into a 1-row data.table
        # Coerce all scalars to character to avoid type clashes, then convert numerics later if you wish
        t2 <- as.data.table(as.list(lapply(node, function(x) if (length(x) == 1) x else NA)))
        t2[, sample_name := sample]
        setcolorder(t2, c("sample_name", setdiff(names(t2), "sample_name")))
        if (is.null(stats)) stats <- t2 else stats <- rbind(stats, t2, fill = TRUE)
    }
    if (missing_files > 0) warning(sprintf("Stats files missing for %s samples.", missing_files))
    return(stats)
}


getResultFilesFromYaml <- function(samples, results_subdir, field) {
    # field should be "mapped_bam" or "bigwig"
    rf <- NULL
    for (sample in samples) {
        yaml_file <- file.path(results_subdir, sample, "stats.yaml")
        if (!file.exists(yaml_file)) next
        y <- yaml::read_yaml(yaml_file)
        pipe_name <- names(y)[1]
        node <- tryCatch(y[[pipe_name]][["sample"]][[sample]], error = function(e) NULL)
        if (is.null(node)) next
        if (!is.null(node[[field]])) {
            t <- data.table(sample_name = sample,
                            object = field,
                            val = as.character(node[[field]]))
            rf <- if (is.null(rf)) t else rbind(rf, t, fill = TRUE)
        }
    }
    return(rf)
}


################################################################################
##### MAIN #####

# Set the project configuration file
pep <- argv$config
# Load the project
prj <- invisible(suppressWarnings(pepr::Project(pep)))
# Convenience
project_name    <- config(prj)$name
pipeline_mode <- config(prj)$pipeline_mode
project_protocol <- unique(invisible(suppressWarnings(pepr::sampleTable(prj)$protocol)))
project_samples <- pepr::sampleTable(prj)$sample_name

sample_table <- data.table(prj@samples)

# Set the output directory
summary_dir <- suppressMessages(file.path(argv$output, "summary"))
# Produce output directory (if needed)
dir.create(summary_dir, showWarnings = FALSE)

# Get project genomes
genomes <- invisible(suppressWarnings(pepr::sampleTable(prj)$genome))
genome <- unique(genomes)

# Set the results subdirectory
if (dir.exists(argv$results)) {
    results_subdir <- suppressMessages(file.path(argv$results))
} else {
    warning(paste0("The project results subdirectory (", argv$results,
            ") does not exist."))
    quit()
}

# ########################
# # Generate stats summary (now from YAML)
stats <- createStatsSummary(project_samples, results_subdir)
if (is.null(stats) || nrow(stats) == 0) quit()
project_stats_file <- file.path(summary_dir, paste0(project_name, '_stats_summary.tsv'))
fwrite(stats, project_stats_file, sep = "\t", col.names = TRUE)


# ########################
# # Generate BAM files summary and igv_session files
# write(paste0("Creating BAM files summary..."), stdout())
# bam_files <- getResultFilesFromYaml(project_samples, results_subdir, "mapped_bam")
# bam_files_tsv <- copy(bam_files)
# # For compatibility with your previous TSV (project-level absolute path column)
# # If you still want to prepend argv$results, keep it; otherwise paths are already absolute in YAML.
# bam_files_tsv$val <- ifelse(startsWith(bam_files_tsv$val, "/"),
#                             bam_files_tsv$val,
#                             file.path(argv$results, bam_files_tsv$val))

# project_bam_files <- file.path(argv$output, paste0(project_name, '_BAM_files.tsv'))
# message(sprintf("Summary (n=%s): %s", length(unique(stats$sample_name)), project_bam_files))
# fwrite(bam_files_tsv, project_bam_files, sep = "\t", col.names = TRUE)

# # IGV session for BAM
# write(paste0("Creating BAM IGV session file..."), stdout())
# project_bam_igv <- file.path(argv$output, paste0(project_name, '_BAM_igv_session.xml'))
# write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>", project_bam_igv)
# write(paste("<Session genome=\"",genome,"\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"All\" version=\"8\">", sep=""),
#       project_bam_igv, append=TRUE)
# write("\t<Resources>", project_bam_igv, append=TRUE)
# for (row in 1:nrow(bam_files)) {
#   sample_name <- bam_files[row, "sample_name"]
#   bam_file <- bam_files[row, "val"]
#   # Optional path rewrite:
#   bam_file <- gsub("/work/project", "X:", bam_file)
#   write(paste("<Resource path=\"", bam_file, "\"/>", sep=""), project_bam_igv, append=TRUE)
# }
# write("\t</Resources>", project_bam_igv, append=TRUE)
# write("</Session>", project_bam_igv, append=TRUE)

# ########################
# # Generate BigWig igv_session file
# write(paste0("Creating BigWig IGV session file..."), stdout())
# bw_files <- getResultFilesFromYaml(project_samples, results_subdir, "bigwig_dedup")
# project_bw_igv <- file.path(argv$output, paste0(project_name, '_BigWig_igv_session.xml'))
# write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>", project_bw_igv)
# write(paste("<Session genome=\"",genome,"\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"All\" version=\"8\">", sep=""),
#       project_bw_igv, append=TRUE)
# write("\t<Resources>", project_bw_igv, append=TRUE)
# for (row in 1:nrow(bw_files)) {
#   bw_file <- bw_files[row, "val"]
#   bw_file <- gsub("/work/project", "X:", bw_file)   # same optional rewrite
#   write(paste("<Resource path=\"", bw_file, "\"/>", sep=""), project_bw_igv, append=TRUE)
# }
# write("\t</Resources>", project_bw_igv, append=TRUE)
# write("</Session>", project_bw_igv, append=TRUE)


#######################################
# Generate FeatureCount classes summary
if (pipeline_mode == "repeats")
{
  write(paste0("Creating feature count classes summary..."), stdout())

  feature_counts_file <- file.path(summary_dir,
                                  paste0(project_name, '_fc_summary.tsv'))
  feature_counts_rds <- file.path(summary_dir,
                                   paste0(project_name, '_fc_summary.rds'))
  for (sample in project_samples) {
    sample_output_folder <- file.path(results_subdir, sample)
    sample_fc_file   <- file.path(sample_output_folder, "feature_counts",
                                  paste(sample,".fc.txt",sep=""))
    t <- fread(sample_fc_file, header=F,
               col.names=c('repeatID', 'length', sample), skip = 2)
    if (exists("fc", inherits = F)) {
      fc <- cbind(fc, t[,3])
    } else {
      fc <- t
    }
  }

  #fwrite(fc, feature_counts_file, sep="\t", col.names=TRUE)

  #generate SummarizedExperiment
  fcm <- as.matrix(fc[,3:ncol(fc)])
  rownames(fcm) <- fc$repeatID
  fc.se <- SummarizedExperiment(assays = list(counts=fcm), colData = sample_table)
  saveRDS(fc.se, file = feature_counts_rds)
}

###################################################
# Generate FeatureCount individual elements summary
if (pipeline_mode == "repeats")
{
  write(paste0("Creating feature counts elements summary..."), stdout())

  feature_counts_id_file <- file.path(summary_dir,
                                  paste0(project_name, '_fc_id_summary.tsv'))
  feature_counts_id_rds <- file.path(summary_dir,
                                   paste0(project_name, '_fc_id_summary.rds'))
  for (sample in project_samples) {
    sample_output_folder <- file.path(results_subdir, sample)
    sample_fc_file   <- file.path(sample_output_folder, "feature_counts",
                                  paste(sample,".fc.id.txt",sep=""))
    t <- fread(sample_fc_file, header=F,
               col.names=c('repeatID', 'chr', 'start', 'end', 'strand', 'length', sample), skip = 2)
    if (exists("fcid", inherits = F)) {
      fcid <- cbind(fcid, t[,7])
    } else {
      fcid <- t
    }
  }

  #generate SummarizedExperiment
  fcidm <- as.matrix(fcid[,7:ncol(fcid)])
  rownames(fcidm) <- fcid$repeatID
  fcid.se <- SummarizedExperiment(assays = list(counts=fcidm), colData = sample_table)
  saveRDS(fcid.se, file = feature_counts_id_rds)
}

###################################################
# Generate IAP coverage summary (for mouse samples)
if (genome =="mm10" & pipeline_mode == "repeats")
{
	write(paste0("Creating IAP coverage summary..."), stdout())

	IAP_coverage_file <- file.path(summary_dir,
                                 paste0(project_name, '_IAP_coverage_summary.tsv'))
  IAP_coverage_rds <- file.path(summary_dir,
                                 paste0(project_name, '_IAP_coverage_summary.rds'))

	for (sample in project_samples) {
	  sample_output_folder <- file.path(results_subdir, sample)
	  sample_IAP_coverage_file   <- file.path(sample_output_folder, "IAP_coverage",
	                                paste(sample,".IAP.norm.coverage.txt",sep=""))
	  cf <- fread(sample_IAP_coverage_file, header=F,
        	     col.names=c(sample))
	  if (exists("coverage", inherits = F)) {
	    coverage <- cbind(coverage, cf)
	  } else {
	    cf$pos <- c(1:nrow(cf))
	    coverage <- cf[,c(2,1)]
	  }
	}
	#fwrite(coverage, file = IAP_coverage_file, sep="\t",
  #      	    col.names=TRUE, row.names=F, quote=F)

#generate SummarizedExperiment
cm <- as.matrix(coverage[,2:ncol(coverage)])
rownames(cm) <- coverage$pos
cov.se <- SummarizedExperiment(assays = list(IAP.coverage=cm), colData = sample_table)
saveRDS(cov.se, file = IAP_coverage_rds)


	# Generate IAP coverage plots
	write(paste0("Creating IAP coverage plots..."), stdout())

	IAP_coverage_merged <- file.path(summary_dir,
        	                       paste0(project_name, '_IAP_coverage_merged.pdf'))
	IAP_coverage_samples <- file.path(summary_dir,
	                               paste0(project_name, '_IAP_coverage_samples.pdf'))
	IAP_coverage_merged_png <- file.path(summary_dir,
        	                         paste0(project_name, '_IAP_coverage_merged.png'))
	IAP_coverage_samples_png <- file.path(summary_dir,
        	                          paste0(project_name, '_IAP_coverage_samples.png'))


	df <- melt(coverage ,  id.vars = 'pos', variable.name = 'samples')

	gm <- ggplot(df, aes(pos,value)) + geom_line(aes(colour = samples))
	ggsave(gm, filename = IAP_coverage_merged)
	ggsave(gm, filename = IAP_coverage_merged_png)

	gs <- ggplot(df, aes(pos,value)) + geom_line() + facet_grid(samples ~ .)
	ggsave(gs, filename = IAP_coverage_samples)
	ggsave(gs, filename = IAP_coverage_samples_png)
}

################################################################################
#Summarize Gene Counts for RNA-seq data (use column 3 with stranded information)
#select RNA-seq samples
RNA_samples <- sample_table[sample_table$protocol == "RNA", "sample_name"]

#if RNA samples were processed
if (nrow(RNA_samples) > 0) {
	write(paste0("Creating gene counts summary..."), stdout())

	for (sample in RNA_samples$sample_name) {
		sample_output_folder <- file.path(results_subdir, sample)
		sample_gc_file   <- file.path(sample_output_folder, paste("aligned_",genome,sep=""),
                                paste(sample,".ReadsPerGene.out.tab",sep=""))
		t <- fread(sample_gc_file, header=F,
				col.names=c('geneID', sample, 'read1', 'read2'), skip = 4)
	  #unstranded (col 2)
	  if (exists("unstranded", inherits = F)) {
		  unstranded <- cbind(unstranded, t[,2])
	  } else {
		  unstranded <- t[,c(1,2)]
	  }
	  #read1 (col 3)
	  colnames(t) <- c('geneID', 'unstranded', sample, 'read2')
	  if (exists("read1", inherits = F)) {
		    read1 <- cbind(read1, t[,3])
		  } else {
		    read1 <- t[,c(1,3)]
	  }
	  #read2 (col 4)
	  colnames(t) <- c('geneID', 'unstranded', 'read1', sample)
	  if (exists("read2", inherits = F)) {
		    read2 <- cbind(read2, t[,4])
		  } else {
		    read2 <- t[,c(1,4)]
	  }
  }
	
  #estimate sense reads which should have more coverage
  if (sum(read1[,2:ncol(read1)]) > sum(read2[,2:ncol(read2)])) {sense <- read1} else {sense <- read2}
	
  #Summarized Experiment
	#unstranded
  unstranded_gene_counts_rds <- file.path(summary_dir,
                                paste0(project_name, '_unstranded_gene_counts_summary.rds'))
	gcm <- as.matrix(unstranded[,2:ncol(unstranded)])
	rownames(gcm) <- unstranded$geneID
	gc.se <- SummarizedExperiment(assays = list(counts=gcm), colData = sample_table[sample_table$sample_name %in% RNA_samples$sample_name,])
	saveRDS(gc.se, file = unstranded_gene_counts_rds)
	#sense reads
	sense_gene_counts_rds <- file.path(summary_dir,
	                                        paste0(project_name, '_sense_gene_counts_summary.rds'))
	gcm <- as.matrix(sense[,2:ncol(sense)])
	rownames(gcm) <- sense$geneID
	gc.se <- SummarizedExperiment(assays = list(counts=gcm), colData = sample_table[sample_table$sample_name %in% RNA_samples$sample_name,])
	saveRDS(gc.se, file = sense_gene_counts_rds)
}
