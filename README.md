# NGS.analysis

NGS.analysis is a pipeline for mapping and quality control of NGS datasets with focus on repeat detection. It is based on the [pypiper](http://pypiper.databio.org/en/latest/) framework and uses [looper](http://looper.databio.org/en/latest/) to submit jobs to a high performance computing cluster. In the ***repeats*** mode mapping and repeat coverage follow the guidelines from [Teissandier, 2019](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1). In ***genes*** mode the pipeline uses bowtie2 (or STAR for RNA-seq) for mapping.

The following experiments (protocols) are currently supported: RNA-seq, ChIP-seq, Cut&Run, Cut&Tag, ATAC-sequencing

The pipeline will produce quality control metrics, BAM and BigWig as well as gene and repeats coverage files for downstream analyses. Basic calculation of differential expression/enrichment can be performed with the NGS.shiny.app. In addition, the pipeline generates lists of the generated BAM files which can be used as sample files for downstream scripts. IGV session files for BAM and BigWig files allow easy visualization.

## Usage

Samples need to be specified in a sample list containing the columns exactly as specified in the example sample file. Don't use spaces or special characters for any entry. The protocol column should be one of CHIP (ChIP-seq), RNA (RNA-seq), CR (Cut&Run), CT (Cut&Tag), ATAC (ATAC-seq).

The project yaml file is used to specify the project name, the pipeline run mode (genes/repeats), the location of the fastq files and additional configuration settings (which usually don't need to be changed).

The pipeline can be started with:
 ```looper run -p slurm project.yaml```

After completion, the summary can be calculated with:
 ```looper runp -p slurm project.yaml```

Finally the reporting summary is generated with:
 ```looper report project.yaml```

## Installation

The pipeline requires installation in a fresh conda environment. The following steps are necessary:

***first use of conda:*** install Anaconda3 from https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
initialize conda base environment:
    conda init

update conda packages:

    conda update -n base -c defaults conda

add conda-forge and bioconda channels:

    conda config --add channels conda-forge
    conda config --add channels bioconda

Create fresh conda environment for ngs pipeline:

    conda env create -n ngs Python=3.9
    conda activate ngs

Install python packages

    conda install -c bioconda samtools
    conda install -c bioconda bedtools
    conda install -c bioconda -c conda-forge bowtie2
    conda install -c bioconda fastqc
    conda install -c bioconda subread
    conda install -c bioconda trimmomatic
    conda install -c bioconda picard
    conda install -c conda-forge -c bioconda deeptools
    conda install -c bioconda star

Install R libraries

open terminal and start R

    install.packages("argparser")
    install.packages("pepr")
    install.packages("reshape2")
    install.packages("ggplot2")
    install.packages("optigrab")

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("SummarizedExperiment")

Install Looper

    python -m pip install looper
    python -m pip install piper

clone NGS.analysis pipeline (in /home folder)

    git clone https://github.com/GunnarSchotta/NGS.analysis

adjust path and configuration settings in:
* NGS.analysis.yaml
* sample_pipeline_interface.yaml
* project_pipeline_interface.yaml
