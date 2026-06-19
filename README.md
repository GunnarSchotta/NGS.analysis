# NGS.analysis

NGS.analysis is a pipeline for mapping and quality control of NGS datasets with focus on repeat detection. It is based on the [pypiper](http://pypiper.databio.org/en/latest/) framework and uses [looper](http://looper.databio.org/en/latest/) to submit jobs to a high performance computing cluster. In the ***repeats*** mode mapping and repeat coverage follow the guidelines from [Teissandier, 2019](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1). In ***genes*** mode the pipeline uses bowtie2 (or STAR for RNA-seq) for mapping.

The following experiments (protocols) are currently supported: RNA-seq, ChIP-seq, Cut&Run, Cut&Tag, ATAC-sequencing

The pipeline will produce quality control metrics, BAM and BigWig as well as gene and repeats coverage files for downstream analyses. Basic calculation of differential expression/enrichment can be performed with the NGS.shiny.app. In addition, the pipeline generates lists of the generated BAM files which can be used as sample files for downstream scripts. IGV session files for BAM and BigWig files allow easy visualization.

## Usage

Samples need to be specified in a sample list containing the columns exactly as specified in the example sample file. Don't use spaces or special characters for any entry. The protocol column should be one of CHIP (ChIP-seq), RNA (RNA-seq), CR (Cut&Run), CT (Cut&Tag), ATAC (ATAC-seq).

The project yaml file is used to specify the project name, the pipeline run mode (genes/repeats), the location of the fastq files and additional configuration settings (which usually don't need to be changed).

The pipeline can be started with:
```
looper run -p slurm
```

After completion, the summary can be calculated with:
```
looper runp -p slurm
```

Finally, generate the custom report from the looper project directory
(the directory that contains `.looper.yaml`):
 ```../generate_report.py```

This writes the main page to:
 ```<results_dir>/report.html```

The script reads `output_dir` from `.looper.yaml`, aggregates directly from the
per-record `stats.yaml` files in the sample and project subdirectories, and
only shows values and objects that were actually reported for the current
pipeline variant.

## Installation

### 1. Install micromamba

Follow the instructions at https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### 2. Clone the pipeline

```bash
git clone https://github.com/GunnarSchotta/NGS.analysis
cd NGS.analysis
```

### 3. Create the environment

```bash
micromamba env create -f environment.yaml
micromamba activate ngs.v2
```

This installs all required tools in one step:
- Alignment: STAR, Bowtie2, Samtools, Bedtools
- Quantification: featureCounts (Subread), RSEM, deepTools
- Quality control: FastQC, Trimmomatic, Picard
- Pipeline framework: looper, piper, pipestat, peppy
- R packages: pepr, data.table, ggplot2, SummarizedExperiment, and others

### 4. Install r-optigrab from CRAN

One R package used for TSS enrichment plots (`plot.TSS.enrichment.R`) is not on
conda-forge and must be installed from CRAN after activating the environment:

```bash
Rscript -e 'install.packages("optigrab", repos="https://cloud.r-project.org")'
```

### 5. Configure the pipeline

Adjust the paths in the two pipeline interface files to match your installation:

- `sample_pipeline_interface.yaml` — path to `NGS.analysis.py`
- `project_pipeline_interface.yaml` — path to `NGS.analysis.collator.py`

Genome indexes and annotation paths are set per-project in `analysis.configuration.yaml`
(see `example/analysis.configuration.yaml`).

### Reproducibility

`requirements.yaml` contains the full pinned environment export used during development.
To reproduce the exact environment:

```bash
micromamba env create -f requirements.yaml
```
