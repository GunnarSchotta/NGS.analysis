library(optigrab)
library(ggplot2)
library(data.table)

#' Plot TSS enrichment
#'
#' This function plots the global TSS enrichment.
#'
#' @param TSSfile TSS enrichment file
#' @param cutoff A TSS enrichment score cutoff value for low quality samples
#' @keywords TSS enrichment
#' @export
#' @examples
#' data("tss")
#' plotTSS(TSSfile = "tss")
#' @export
plotTSS <- function(TSSfile, cutoff=6) {
    if (length(TSSfile) == 1) {
        write(paste0("\nGenerating TSS plot with ", TSSfile), stdout())
    } else {
        if (length(TSSfile) == 2) {
            write(paste0("\nGenerating TSS plot with ",
                         paste(TSSfile, collapse=" and ")),
                  stdout())
        } else {
            write(paste0("\nNot sure how to merge the following: ",
                         paste(TSSfile, collapse=", ")),
                  stdout())
            write(paste0("Did you mean to pass more than 2 files?"), stdout())
            quit(save = "no", status = 1, runLast = FALSE)
        }
    }

    t1 <- theme_classic(base_size=14) +
            theme(plot.background  = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border     = element_rect(colour = "black",
                                                  fill=NA, size=0.5),
                  panel.background = element_blank(),
                  axis.line    = element_blank(),
                  legend.position="none",
                  aspect.ratio = 1,
                  axis.ticks.length = unit(2, "mm"))

    iMat <- data.table(V1 = numeric())
    if (length(TSSfile) == 1) {
        if (exists(TSSfile)) {
            iMat <- data.table(get(TSSfile))
        } else {
            iMat <- fread(TSSfile)
        }
    } else if (length(TSSfile) == 2) {
        for (i in 1:length(TSSfile)) {
            if (exists(TSSfile[i])) {
                if (i == 1) {
                    iMat <- data.table(get(TSSfile[i]))
                } else {
                    iMat <- list(iMat, data.table(get(TSSfile[i])))
                }
            } else {
                if (i == 1) {
                    iMat <- fread(TSSfile[i])
                } else {
                    iMat <- list(iMat, fread(TSSfile[i]))
                }
            }
        }
    } else {
        write(paste0("\nNot sure how to merge the following: ",
                     paste(TSSfile, collapse=", ")),
              stdout())
        write(paste0("Did you mean to pass more than 2 files?"), stdout())
        quit(save = "no", status = 1, runLast = FALSE)
    }

    if (length(TSSfile) == 1) {
        plusMinus <- iMat
    } else {
        plus      <- iMat[[1]]
        minus     <- iMat[[2]]
    }

    if (exists("plusMinus")) {
        val      <- 0.05*nrow(plusMinus)
        #normTSS  <- (plusMinus / mean(plusMinus[c(1:val,
        #            (nrow(plusMinus)-val):nrow(plusMinus)), V1]))
        normTSS           <- plusMinus / mean(plusMinus[c(1:val), V1])
        colnames(normTSS) <- c("score")
        peakPos  <- which.max(normTSS$score)
        # check for true peak
        if ((normTSS$score[peakPos]/normTSS$score[peakPos-1]) > 1.5 &
            (normTSS$score[peakPos]/normTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- normTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile, "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
        }
    } else {
        val      <- 0.05*nrow(plus)
        #normTSS  <- (plus / mean(plus[c(1:val,
        #            (nrow(plus)-val):nrow(plus)), V1]))
        normTSS           <- plus / mean(plus[c(1:val), V1])
        colnames(normTSS) <- c("score")
        peakPos  <- which.max(normTSS$score)
        # check for true peak
        if ((normTSS$score[peakPos]/normTSS$score[peakPos-1]) > 1.5 &
            (normTSS$score[peakPos]/normTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- normTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }
        TSSscore <- round(mean(normTSS[(max(0, peakPos-50)):(min(nrow(normTSS),
                                       peakPos+50)), score]),1)
        if (is.nan(TSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[1], "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
        }
    }

    lineColor <- "red2"
    if (TSSscore > cutoff)
    {
        lineColor <- "springgreen4"
    }

    name        <- basename(tools::file_path_sans_ext(TSSfile[1]))
    numFields   <- 2
    for(j in 1:numFields) name <- gsub("_[^_]*$", "", name)
    sample_name <- paste(dirname(TSSfile[1]), name, sep="/")

    pre <- ggplot(normTSS, aes(x=(as.numeric(rownames(normTSS))-
                                 (nrow(normTSS)/2)),
                               y=score, group=1, colour="black")) +
        # geom_hline(yintercept = 6, linetype = 2,
        #            color = "grey", size = 0.25) +
        geom_smooth(method="loess", span=0.02,
                    se=FALSE, colour=lineColor) +
        labs(x = "Distance from TSS (bp)", y = "TSS enrichment score")
    y_max <- TSSscore
    p <- pre + t1 +
         scale_x_continuous(expand=c(0,0)) +
         scale_y_continuous(expand=c(0,0)) +
         coord_cartesian(xlim=c(-2300, 2300), ylim=c(0, 1.1*y_max))
    if (exists("minus")) {
        val      <- 0.025*nrow(minus)
        # normTSS  <- (minus / mean(minus[c(1:val,
        #             (nrow(minus)-val):nrow(minus)), V1]))
        minusNormTSS           <- minus / mean(minus[c(1:val), V1])
        colnames(minusNormTSS) <- c("score")
        peakPos       <- which.max(minusNormTSS$score)
        # check for true peak
        if ((minusNormTSS$score[peakPos]/minusNormTSS$score[peakPos-1]) > 1.5 &
            (minusNormTSS$score[peakPos]/minusNormTSS$score[peakPos+1]) > 1.5) {
            tmpTSS  <- minusNormTSS$score[-peakPos]
            peakPos <- which.max(tmpTSS) + 1
        }

        minusTSSscore <- round(
            mean(minusNormTSS[(max(0, peakPos-50)):(min(nrow(minusNormTSS),
                               peakPos+50)), score]),1)
        if (is.nan(minusTSSscore)) {
            message(paste0("\nNaN produced.  Check ", TSSfile[2], "\n"))
            quit(save = "no", status = 1, runLast = FALSE)
        }
        p <- p + geom_smooth(data=minusNormTSS,
                             aes(x=(as.numeric(rownames(minusNormTSS))-
                                   (nrow(minusNormTSS)/2)),
                                 y=score, group=1, colour="black"),
                             method="loess", span=0.02,
                             se=FALSE, colour="blue") +
                annotate("rect", xmin=1200, xmax=2300, ymin=0.9*y_max,
                         ymax=1.1*y_max, fill="gray95") +
                annotate("text", x=1750, y=1.05*y_max, label="TSS Score",
                         fontface = 1, hjust=0.5) +
                annotate("text", x=1500, y=y_max, label="+", fontface = 2,
                          hjust=0.5, color=lineColor) +
                annotate("text", x=1500, y=0.95*y_max, label=TSSscore,
                         fontface = 2,  hjust=0.5, color=lineColor) +
                annotate("text", x=2000, y=y_max, label="-",
                         fontface = 2,  hjust=0.5, color="blue") +
                annotate("text", x=2000, y=0.95*y_max, label=minusTSSscore,
                         fontface = 2,  hjust=0.5, color="blue")
    } else {
        p <- p + annotate("rect", xmin=1200, xmax=2300, ymin=0.9*y_max,
                          ymax=1.1*y_max, fill="gray95") +
                 annotate("text", x=1750, y=1.05*y_max, label="TSS Score",
                          fontface = 1, hjust=0.5) +
                 annotate("text", x=1750, y=0.95*y_max, label=TSSscore,
                          fontface = 2, hjust=0.5)
    }

    return(p)
}


TSS_CUTOFF <- 6
# get arguments
a   <- opt_get_args()
# Determine the number of input files
inArgs <- 0
p      <- 1
val    <- a[p]
while (!(val %in% c("-i", "--input"))) {
    p   <- p + 1
    val <- a[p]
}
while (p < length(opt_get_args())) {
    p   <- p + 1
    val <- a[p]
    if (!(val %in% c("-i", "--input"))) {
        inArgs <- inArgs + 1
    }
}


#' Derive the sample name from input file and return with full path
#'
#' @param path A path to a file for which you wish to extract the sample name
#' @param num_fields An integer representing the number of fields to strip
#' @param delim A delimiter for the fields splitting a path or string
#'
#' @export
sampleName <- function(path, num_fields=2, delim='_') {
    name <- basename(tools::file_path_sans_ext(path))
    if(num_fields == 0) {return(name)}
    for(n in 1:num_fields) name <- gsub(paste0(delim, "[^", delim, "]*$"), "", name)
    return(paste(dirname(path), name, sep="/"))
}

#plot TSS based on TSS enrichment file

TSSfile     <- opt_get(name = c("input", "i"), required=TRUE, n=inArgs,
                       description="TSS enrichment file.")
print(TSSfile)
p           <- plotTSS(TSSfile = TSSfile)

sample_name <- sampleName(TSSfile[1])

png(filename = paste0(sample_name, "_TSS_enrichment.png"),
    width = 275, height = 275)
suppressWarnings(print(p))
invisible(dev.off())

pdf(file = paste0(sample_name, "_TSS_enrichment.pdf"),
    width = 4, height = 4, useDingbats=F)
suppressWarnings(print(p))
invisible(dev.off())

if (exists("p")) {
    write("TSS enrichment plot completed!\n", stdout())
} else {
    write("Unable to produce TSS enrichment plot!\n", stdout())
}
