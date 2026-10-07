# NGS.analysis (v3)

NGS.analysis maps and QCs ChIP-seq, ATAC-seq, CUT&RUN (CR), CUT&Tag (CT) and RNA-seq data, with a focus on repeats.
It is built on [pypiper](http://pypiper.databio.org/en/latest/) and submitted to SLURM with
[looper](http://looper.databio.org/en/latest/).

**Changes from v2:** see `CHANGELOG.md`.

**Modes**
- **`genes`:**
  - chromatin: bowtie2 → ENCODE-style filtered BAM → CPM bigwig, peaks + FRiP, QC;
  - RNA: STAR + RSEM.
- **`repeats`:** mapping and repeat quantification follow [Teissandier et al. 2019](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1). This is unchanged from v2; open questions are in `docs/repeat_strategy_open_questions.md`.
  - STAR reports one random locus per multimapper.
  - featureCounts counts at family level (all reads) and element level (unique reads).
  - IAP gag coverage is computed for mm10.

**Protocols** (`protocol` column): CHIP, ATAC, CR, CT, RNA.

## What the pipeline does

| Step | genes mode | repeats mode |
|---|---|---|
| Trimming (all modes) | Trimmomatic `ILLUMINACLIP:<Nextera or TruSeq>:2:30:10:1:true MINLEN:30` (PE: adapter overhangs clipped down to 1 nt, both reads of short fragments kept); SE: `:2:30:7` | same |
| Alignment, chromatin | `bowtie2 --very-sensitive -X 2000 --dovetail` | STAR, one random locus per multimapper (v2 settings) |
| Alignment, RNA | STAR (default multimapping) + RSEM | STAR multimapper settings + RSEM |
| Duplicates | Picard MarkDuplicates | Picard MarkDuplicates (marked only) |
| Analysis BAM | `<s>.filt.bam`: proper pairs, duplicates **removed**, MAPQ ≥ 30, canonical chromosomes, no chrM, blacklist removed, orphan mates removed | `<s>.dedup.unique.bam`: MAPQ 255, duplicates **marked, not removed** (v2) |
| Bigwig | `<s>.filt.bw` (CPM, fragment-extended); RNA: `<s>.bw` | `<s>.dedup.unique.bw` (RPKM, v2) |
| CR/CT | `.filt.nuc/.subnuc` BAMs + bigwigs (fragments ≥ / < 120 bp) | `.dedup.unique.nuc/.subnuc` (v2) |
| QC | insert size, library complexity (NRF/PBC), TSS score (ATAC), mito %, duplication | same, on the repeats BAM |
| Peaks | MACS3, narrow/broad (`peak_mode`, `auto` from `target`), blacklist-filtered, FRiP | same, on the unique BAM |
| Repeats | – | featureCounts family / element counts, IAP coverage (mm10) |

All statistics are in **fragments** (PE pairs, SE reads), and percentages are on a 0–100 scale. Every key is
defined in `docs/stats_units.md`.

## Usage

1. **Project folder.** Create `sample.table.csv`, `analysis.configuration.yaml` and `.looper.yaml` (see `example/`).
   - **Sample table:** `sample_name, ..., protocol, organism, fastq1, fastq2, read_type, read1, read2`.
   - **Optional columns:**
     - `target`: antibody; used by `peak_mode: auto`. Broad for H3K27me3/H3K9me3/H3K36me3/H3K79me2/H4K20me3; no peaks for none/input/IgG.
     - `peak_mode`: narrow, broad, none or auto.
     - `control`: sample_name of the input/IgG, used by `NGS.peaks`.
   - **Sample names:** don't use spaces or special characters, and don't let one sample name be a prefix of another at a `_` boundary (`A` and `A_rep2`). pipestat's status-flag lookup cannot tell them apart.
   - **Genome:** set by PEP `imply` from `organism` (`genome: mm10 | hg38 | hs1 | mhaESC`). Indices and annotation come from `genomes/<genome>.yaml`; `genomes/README.md` lists what is available.
2. **Preflight check** in the project folder:
   ```
   /store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/python3 <v3>/check_project.py
   ```
3. **Run:**
   ```
   export PATH=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin:$PATH
   looper run -p slurm            # per sample
   looper runp -p slurm           # project summary, after all samples have finished (add --ignore-flags to re-collate)
   python <v3>/generate_report.py # results/report.html
   ```
4. **Peaks against a control** (ChIP vs input/IgG), after the main run. Use a copy of `.looper.yaml` whose `pipeline_interfaces` lists only `<v3>/peaks_pipeline_interface.yaml`:
   ```
   looper run -p slurm -c .looper.peaks.yaml
   ```

v3 refuses to write into an output folder that contains v2 results, because pypiper would silently reuse old files.
Use a new `output_dir`.

## Project summary (`results/project/summary/`)

- `<project>_stats_summary.tsv`: all statistics; file objects as absolute paths.
- `<project>_files.tsv`, `<project>_BAM_files.tsv`: output files per sample.
- repeats mode: `_fc_summary.rds`, `_fc_id_summary.rds` (SummarizedExperiment) and IAP coverage (mm10).
- RNA: `_unstranded_` / `_sense_gene_counts_summary.rds`, `_tpm_summary.tsv/.rds`.
- `results/<project>_BigWig_igv_session.xml`.
- The genome is appended to file names when a project contains more than one genome.

## Installation

The pipeline runs in the `ngs.v3` micromamba environment (`/store24/project24/becgsc_001/micromamba/envs/ngs.v3`).
- `env/build_ngs.v3.sh` rebuilds it: an exact copy of `ngs.v2` (conda packages, pip packages, R packages) plus MACS3 3.0.5 and gffread.
- The pinned definitions are `env/ngs.v3.explicit.txt` and `env/ngs.v3.pip.txt`.
- The interfaces call the environment's python by absolute path and find scripts via `{looper.piface_dir}`, so no paths need to be edited.

## Validation and benchmark

- `tests/`: T1 smoke-test data and run scripts (v2 vs v3), `test_schema_keys.py`, `compare_T1.py`.
- `benchmark/`: ART simulation (per TE family mapping % and true-positive rate) and an adapter read-through simulation (fragment recovery by length).
- Results are in `/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/` and `.../NGS.analysis.v3.benchmark/`.
