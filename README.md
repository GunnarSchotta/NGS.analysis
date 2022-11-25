# NGS.repeats.analysis

NGS.repeats.analysis is a pipeline for mapping and quality control of NGS datasets with focus on repeat detection. Mapping and repeat coverage follow the guidelines from [Teissandier, 2019](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1).

The following experiments (protocols) are currently supported: RNA-seq, ChIP-seq, Cut&Run, Cut&Tag, ATAC-sequencing

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

clone NGS.repeats.pipeline (in /home folder)

    git clone https://github.com/gscarnuntum/NGS.repeats.analysis

adjust path and configuration settings in:
* NGS.repeats.analysis.yaml
* sample_pipeline_interface.yaml
* project_pipeline_interface.yaml
