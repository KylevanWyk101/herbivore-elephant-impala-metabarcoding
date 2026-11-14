# Dual-Locus Metabarcoding: Elephant vs. Impala Dry-Season Diet

Code behind: van Wyk, Botha, Barnard & Siebert (2026), *Ecology and Evolution*,
"Dual-Locus Metabarcoding Reveals Dry-Season Dietary Differentiation Between Elephants
and Impalas in a South African Savanna." https://doi.org/10.1002/ece3.74487

Raw sequencing data: NCBI BioProject [PRJNA1465205](http://www.ncbi.nlm.nih.gov/bioproject/1465205)
(not included in this repo — see below).

## Pipeline

1. **Sequence processing** — primer trimming (Cutadapt), QC (FastQC/MultiQC), DADA2
   denoising via QIIME2 (`01_sequence_processing/`)
2. **ASV curation** — LULU post-clustering, frequency/abundance filtering (`02_asv_curation/`)
3. **Reference database** — rbcL/trnL libraries built on Botha et al. (2023), expanded with
   58 GenBank accessions and formatted as BLAST databases (`03_reference_database/`)
4. **Taxonomic assignment** — BLASTn against the reference libraries, contamination removal,
   expert curation, final diet-matrix construction (`04_taxonomic_assignment/`)
5. **Diversity statistics** — Shannon index, NMDS/PERMANOVA/PERMDISP, Pianka overlap
   (`05_diversity_statistics/`)
6. **Differential abundance** — DESeq2 by genus (`06_differential_abundance/`)
7. **Visualization** — RRA/FOO bar charts, Sankey diagrams (`07_composition_visualization/`)

## What's not in this repo

Raw FASTQ and QIIME2 intermediate artifacts (~14GB) are excluded — see the BioProject link
above. Primer trimming and QC (Cutadapt, FastQC/MultiQC) were run as one-off commands and
aren't preserved as standalone scripts. The full BLAST-querying/ecologist-review workflow
used to build the reference database (genus/family lookup, confidence scoring) exists but
isn't included here — only the finished, BLAST-formatted databases are. The pipeline order
above reflects what was actually run.

## Dependencies

**Python**: pandas, numpy, matplotlib, seaborn, scipy, scikit-bio, scikit-learn, statsmodels,
biopython, Pillow, plotly. See `environment/requirements.txt`.

**R**: lulu, DESeq2, ggplot2, ggrepel, dplyr, pgirmess, vegan. See `environment/r_packages.txt`.

Sequence processing additionally requires: Cutadapt, FastQC, MultiQC, QIIME2
(`qiime2-amplicon-2024.10`), BLAST+ (`makeblastdb`/`blastn`).

## Reproducibility note

Several R scripts (LULU curation, Pianka index, DESeq2) have hardcoded `setwd()` paths from
the original development machine. Change these to relative paths (or set your working
directory to the script's own folder) before re-running.

