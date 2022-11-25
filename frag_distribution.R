#! /usr/bin/env Rscript
#
# frag_distribution.R
#
# usage: Rscript /path/to/Rscript/frag_distribution.R 
#
# Depends: R (>= 3.5.1)
# Imports: PEPATACr, argparser
###############################################################################
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

# Create a parser
p <- arg_parser("Produce Fragment distribution plot")

# Add command line arguments
p <- add_argument(p, arg="sample", short="-s",
                  help="Sample name")
p <- add_argument(p, arg="fraglen", short="-f",
                  help="Fragment length file")
p <- add_argument(p, arg="output", short="-o",
                  help="Project parent output directory path")

# Parse the command line arguments
argv <- parse_args(p)

fraglen.file <- argv$fraglen
od <- argv$output
sample <- argv$sample

of.pdf <- paste(od,"/",sample,".fraglen.pdf", sep="")
of.png <- paste(od,"/",sample,".fraglen.png", sep="")

fl <- read.table(fraglen.file)
fl$V1 <- abs(fl$V1)

pdf(file=of.pdf)
hist(fl$V1, nclass=50, xlab = "insert size", main = sample)
dev.off()

png(file=of.png)
hist(fl$V1, nclass=50, xlab = "insert size", main = sample)
dev.off()
