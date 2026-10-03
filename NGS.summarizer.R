#! /usr/bin/env Rscript
#
# NGS.summarizer.R (NGS.analysis v3)
#
# Project-level summary files, called by NGS.analysis.collator.py:
#   Rscript NGS.summarizer.R <project_config.yaml> <results/project> <results/samples>
# Outputs (results/project/summary/):
#   {project}_stats_summary.tsv        all per-sample stats; file objects as their absolute path
#   {project}_files.tsv                sample, object, genome, absolute path (BAMs, bigwigs, peaks, ...)
#   {project}_BAM_files.tsv            sample, genome, BAM_mapped / BAM_filtered / BAM_dedup_unique
#   {project}[_{genome}]_fc_summary.rds / _fc_id_summary.rds        repeats mode, one per genome
#   {project}_IAP_coverage_summary.rds (+ plots)                     repeats mode, mm10 samples
#   {project}[_{genome}]_unstranded/_sense_gene_counts_summary.rds  RNA samples
#   {project}[_{genome}]_tpm_summary.tsv/.rds                       RNA samples with RSEM
#   results/{project}[_{genome}]_BigWig_igv_session.xml
# The genome suffix is only added when a project contains more than one genome (v2 file names otherwise).
# IGV path mapping: project config key `igv_path_map` (named list from -> to); default /store24/project24/ -> X:/
###############################################################################
suppressWarnings(suppressPackageStartupMessages({
    library(argparser); library(pepr); library(data.table); library(reshape2)
    library(ggplot2); library(yaml); library(SummarizedExperiment)
}))

p <- arg_parser("Produce Summary Reports, Files, and Plots")
p <- add_argument(p, arg="config", short="-c", help="project_config.yaml")
p <- add_argument(p, arg="output", short="-o", help="Project parent output directory path")
p <- add_argument(p, arg="results", short="-r", help="Project results output subdirectory path")
p <- add_argument(p, arg="--new-start", short="-N", flag=TRUE, help="(compatibility, unused)")
argv <- parse_args(p)

sampleNode <- function(sample, results_subdir) {
    yaml_file <- file.path(results_subdir, sample, "stats.yaml")
    if (!file.exists(yaml_file)) return(NULL)
    y <- yaml::read_yaml(yaml_file)
    tryCatch(y[[names(y)[1]]][["sample"]][[sample]], error = function(e) NULL)
}
asScalar <- function(x) {
    if (is.list(x) && !is.null(x[["path"]])) return(as.character(x[["path"]]))
    if (length(x) == 1) return(x)
    NA
}

prj <- invisible(suppressWarnings(pepr::Project(argv$config)))
cfg <- config(prj)
project_name  <- cfg$name
pipeline_mode <- cfg$pipeline_mode
sample_table  <- data.table(prj@samples)
project_samples <- sample_table$sample_name
genome_of <- setNames(as.character(sample_table$genome), project_samples)
genomes_used <- unique(genome_of)
multi_genome <- length(genomes_used) > 1
fname <- function(g, suffix) paste0(project_name, if (multi_genome) paste0("_", g) else "", suffix)
path_map <- if (!is.null(cfg$igv_path_map)) cfg$igv_path_map else list("/store24/project24/" = "X:/")
mapPath <- function(x) { for (k in names(path_map)) x <- gsub(k, path_map[[k]], x, fixed = TRUE); x }

summary_dir <- file.path(argv$output, "summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
results_subdir <- argv$results
if (!dir.exists(results_subdir)) { warning("results subdirectory missing: ", results_subdir); quit() }

################################################################################
# stats summary + file table
write("Creating stats summary ...", stdout())
stats <- NULL; files <- NULL
for (s in project_samples) {
    node <- sampleNode(s, results_subdir)
    if (is.null(node)) { warning("No stats for sample ", s); next }
    row <- as.data.table(lapply(node, asScalar))
    row[, sample_name := s]
    stats <- if (is.null(stats)) row else rbind(stats, row, fill = TRUE)
    for (k in names(node)) if (is.list(node[[k]]) && !is.null(node[[k]][["path"]])) {
        files <- rbind(files, data.table(sample_name = s, object = k, genome = genome_of[[s]],
                                         path = as.character(node[[k]][["path"]])))
    }
}
if (is.null(stats) || nrow(stats) == 0) quit()
setcolorder(stats, c("sample_name", setdiff(names(stats), "sample_name")))
fwrite(stats, file.path(summary_dir, paste0(project_name, "_stats_summary.tsv")), sep = "\t")
if (!is.null(files)) {
    fwrite(files, file.path(summary_dir, paste0(project_name, "_files.tsv")), sep = "\t")
    bams <- files[object %in% c("BAM_mapped", "BAM_filtered", "BAM_dedup_unique")]
    fwrite(dcast(bams, sample_name + genome ~ object, value.var = "path"),
           file.path(summary_dir, paste0(project_name, "_BAM_files.tsv")), sep = "\t")
    # IGV session(s), one per genome
    for (g in genomes_used) {
        bw <- files[genome == g & object %in% c("BigWig", "BigWig_dedup", "BigWig_filtered")]
        if (nrow(bw) == 0) next
        igv_xml <- file.path(dirname(argv$output), fname(g, "_BigWig_igv_session.xml"))
        writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
                     paste0('<Session genome="', g, '" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">'),
                     "  <Resources>", paste0('    <Resource path="', mapPath(bw$path), '"/>'),
                     "  </Resources>", "</Session>"), igv_xml)
    }
}

samplesOf <- function(g, subset = project_samples) subset[genome_of[subset] == g]
seFrom <- function(m, samples) SummarizedExperiment(assays = list(counts = m),
                                                    colData = sample_table[sample_table$sample_name %in% samples, ])

################################################################################
# repeats mode: featureCounts family / element summaries, IAP coverage
if (pipeline_mode == "repeats") {
    for (g in genomes_used) {
        ss <- samplesOf(g)
        for (lvl in c("fc", "fc_id")) {
            write(paste0("Creating ", lvl, " summary (", g, ") ..."), stdout())
            mat <- NULL; ok <- c()
            for (s in ss) {
                f <- file.path(results_subdir, s, "feature_counts", paste0(s, if (lvl == "fc") ".fc.txt" else ".fc.id.txt"))
                if (!file.exists(f)) { warning("missing ", f); next }
                if (lvl == "fc") {
                    t <- fread(f, header = FALSE, col.names = c("repeatID", "length", s), skip = 2)
                } else {
                    t <- fread(f, header = FALSE, col.names = c("repeatID", "chr", "start", "end", "strand", "length", s), skip = 2)
                }
                mat <- if (is.null(mat)) t else cbind(mat, t[, ncol(t), with = FALSE])
                ok <- c(ok, s)
            }
            if (is.null(mat)) next
            first <- if (lvl == "fc") 3 else 7
            m <- as.matrix(mat[, first:ncol(mat)]); rownames(m) <- mat$repeatID
            saveRDS(seFrom(m, ok), file.path(summary_dir, fname(g, paste0("_", lvl, "_summary.rds"))))
        }
    }
    ss <- samplesOf("mm10")
    cov <- NULL
    for (s in ss) {
        f <- file.path(results_subdir, s, "IAP_coverage", paste0(s, ".IAP.norm.coverage.txt"))
        if (!file.exists(f)) next
        cf <- fread(f, header = FALSE, col.names = s)
        cov <- if (is.null(cov)) data.table(pos = seq_len(nrow(cf)), cf) else cbind(cov, cf)
    }
    if (!is.null(cov)) {
        write("Creating IAP coverage summary ...", stdout())
        cm <- as.matrix(cov[, 2:ncol(cov)]); rownames(cm) <- cov$pos
        cov.se <- SummarizedExperiment(assays = list(IAP.coverage = cm),
                                       colData = sample_table[sample_table$sample_name %in% colnames(cm), ])
        saveRDS(cov.se, file.path(summary_dir, paste0(project_name, "_IAP_coverage_summary.rds")))
        df <- melt(cov, id.vars = "pos", variable.name = "samples")
        gm <- ggplot(df, aes(pos, value)) + geom_line(aes(colour = samples))
        gs <- ggplot(df, aes(pos, value)) + geom_line() + facet_grid(samples ~ .)
        for (ext in c("pdf", "png")) {
            ggsave(gm, filename = file.path(summary_dir, paste0(project_name, "_IAP_coverage_merged.", ext)))
            ggsave(gs, filename = file.path(summary_dir, paste0(project_name, "_IAP_coverage_samples.", ext)))
        }
    }
}

################################################################################
# RNA-seq: STAR gene counts (unstranded / sense) and RSEM TPM, per genome
rna_all <- sample_table[sample_table$protocol == "RNA", ]$sample_name
for (g in genomes_used) {
    ss <- samplesOf(g, rna_all)
    if (length(ss) == 0) next
    write(paste0("Creating gene count summaries (", g, ") ..."), stdout())
    un <- NULL; r1 <- NULL; r2 <- NULL; ok <- c()
    for (s in ss) {
        f <- file.path(results_subdir, s, paste0("aligned_", g), paste0(s, ".ReadsPerGene.out.tab"))
        if (!file.exists(f)) { warning("missing ", f); next }
        t <- fread(f, header = FALSE, col.names = c("geneID", "u", "s1", "s2"), skip = 4)
        un <- if (is.null(un)) t[, .(geneID, u)] else cbind(un, t[, .(u)])
        r1 <- if (is.null(r1)) t[, .(geneID, s1)] else cbind(r1, t[, .(s1)])
        r2 <- if (is.null(r2)) t[, .(geneID, s2)] else cbind(r2, t[, .(s2)])
        ok <- c(ok, s)
    }
    if (length(ok)) {
        for (x in list(un, r1, r2)) setnames(x, c("geneID", ok))
        sense <- if (sum(r1[, -1]) > sum(r2[, -1])) r1 else r2   # project-wide strand guess, as in v2
        for (pair in list(list(un, "_unstranded_gene_counts_summary.rds"), list(sense, "_sense_gene_counts_summary.rds"))) {
            m <- as.matrix(pair[[1]][, -1]); rownames(m) <- pair[[1]]$geneID
            saveRDS(seFrom(m, ok), file.path(summary_dir, fname(g, pair[[2]])))
        }
    }
    tpm <- NULL
    for (s in ss) {
        f <- file.path(results_subdir, s, paste0("aligned_", g), paste0(s, ".rsem.genes.results"))
        if (!file.exists(f)) next
        t <- fread(f, header = TRUE, sep = "\t", select = c("gene_id", "TPM")); setnames(t, "TPM", s)
        tpm <- if (is.null(tpm)) t else merge(tpm, t, by = "gene_id", all = TRUE)
    }
    if (!is.null(tpm)) {
        fwrite(tpm, file.path(summary_dir, fname(g, "_tpm_summary.tsv")), sep = "\t")
        m <- as.matrix(tpm[, -1]); rownames(m) <- tpm$gene_id
        saveRDS(SummarizedExperiment(assays = list(TPM = m),
                                     colData = sample_table[sample_table$sample_name %in% colnames(m), ]),
                file.path(summary_dir, fname(g, "_tpm_summary.rds")))
    }
}
write("Summary done.", stdout())
