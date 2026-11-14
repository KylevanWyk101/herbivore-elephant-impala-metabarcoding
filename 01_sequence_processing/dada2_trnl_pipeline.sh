#!/bin/bash

# Function to display loading message
loading_message() {
    echo -e "\n$1"
    echo "------------------------------------"
    sleep 1
}

# Check if QIIME 2 environment is activated
if [[ -z "$CONDA_DEFAULT_ENV" || "$CONDA_DEFAULT_ENV" != "qiime2-amplicon-2024.10" ]]; then
    echo "Error: QIIME 2 environment (qiime2-amplicon-2024.10) not activated."
    echo "Please activate it with: conda activate qiime2-amplicon-2024.10"
    exit 1
fi
loading_message "QIIME 2 environment verified for trnl analysis"

# Check if manifest file exists
if [[ ! -f "manifest_trnl.tsv" ]]; then
    echo "Error: Manifest file 'manifest_trnl.tsv' not found in current directory."
    exit 1
fi
loading_message "Manifest file verified in current working directory"

# Create output directory
mkdir -p trnl_output
loading_message "Created output directory for trnl results"

# Import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_trnl.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path trnl_output/imported_trnl.qza
loading_message "Imported trnl sequences successfully"

# Visualize imported data
qiime demux summarize \
  --i-data trnl_output/imported_trnl.qza \
  --o-visualization trnl_output/imported_trnl.qzv
loading_message "Generated visualization of imported trnl sequences"

# Denoise with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs trnl_output/imported_trnl.qza \
  --p-trunc-len 0 \
  --p-trunc-q 20 \
  --p-max-ee 2 \
  --p-chimera-method pooled \
  --p-n-threads 4 \
  --o-table trnl_output/trnl_filtered_table.qza \
  --o-representative-sequences trnl_output/trnl_filtered_rep_seqs.qza \
  --o-denoising-stats trnl_output/trnl_stats.qza
loading_message "Processed trnl sequences with DADA2 denoising"

# Export denoising stats to TSV (needed for read count summary)
qiime tools export \
  --input-path trnl_output/trnl_stats.qza \
  --output-path trnl_output/trnl_exported_stats

# Summarize read counts from denoising stats
read_summary=$(awk -F'\t' 'NR==2 {next} NR>1 {input+=$2; filtered+=$3; non_chimeric+=$7} END {printf "Denoising Summary:\n - Total Input Reads: %d\n - Total Filtered Reads: %d\n - Total Non-Chimeric Reads: %d\n", input, filtered, non_chimeric}' E_trnl_output/E_trnl_exported_stats/stats.tsv)
loading_message "$read_summary"

# Visualize denoising stats
qiime metadata tabulate \
  --m-input-file trnl_output/trnl_stats.qza \
  --o-visualization trnl_output/trnl_stats.qzv
loading_message "Generated visualization of trnl denoising statistics"

# Summarize feature table
qiime feature-table summarize \
  --i-table trnl_output/trnl_filtered_table.qza \
  --o-visualization trnl_output/trnl_filtered_table.qzv
loading_message "Summarized trnl feature table"

# Export feature table
qiime tools export \
  --input-path trnl_output/trnl_filtered_table.qza \
  --output-path trnl_output/trnl_exported_table
loading_message "Exported trnl feature table"

# Convert to TSV
biom convert \
  -i trnl_output/trnl_exported_table/feature-table.biom \
  -o trnl_output/trnl_exported_table/trnl_feature_table.tsv \
  --to-tsv
loading_message "Converted trnl feature table to TSV format"

# Export representative sequences
qiime tools export \
  --input-path trnl_output/trnl_filtered_rep_seqs.qza \
  --output-path trnl_output/trnl_fasta_exported
loading_message "Exported trnl representative sequences"

# Rename FASTA
mv trnl_output/trnl_fasta_exported/dna-sequences.fasta trnl_output/trnl_fasta_exported/dada2_trnl.fasta
loading_message "Renamed E_trnl FASTA file to dada2_trnl.fasta"

echo -e "\ntrnl pipeline completed successfully"
