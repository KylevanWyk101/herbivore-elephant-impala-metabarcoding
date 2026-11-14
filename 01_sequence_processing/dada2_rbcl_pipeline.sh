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
loading_message "QIIME 2 environment verified for rbcl analysis"

# Check if manifest file exists
if [[ ! -f "manifest_rbcl.tsv" ]]; then
    echo "Error: Manifest file 'manifest_rbcl.tsv' not found in current directory."
    exit 1
fi
loading_message " Manifest file verified in current working directory"

# Create output directory
mkdir -p rbcl_output
loading_message "Created output directory for rbcl results"

# Import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_rbcl.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path rbcl_output/imported_rbcl.qza
loading_message "Imported rbcl sequences successfully"

# Visualize imported data
qiime demux summarize \
  --i-data rbcl_output/imported_rbcl.qza \
  --o-visualization rbcl_output/imported_rbcl.qzv
loading_message "Generated visualization of imported rbcl sequences"

# Denoise with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs rbcl_output/imported_rbcl.qza \
  --p-trunc-len 0 \
  --p-trunc-q 20 \
  --p-max-ee 2 \
  --p-chimera-method pooled \
  --p-n-threads 4 \
  --o-table rbcl_output/rbcl_filtered_table.qza \
  --o-representative-sequences rbcl_output/rbcl_filtered_rep_seqs.qza \
  --o-denoising-stats rbcl_output/rbcl_stats.qza
loading_message "Processed rbcl sequences with DADA2 denoising"

# Export denoising stats to TSV (needed for read count summary)
qiime tools export \
  --input-path rbcl_output/rbcl_stats.qza \
  --output-path rbcl_output/rbcl_exported_stats

# Summarize read counts from denoising stats
read_summary=$(awk -F'\t' 'NR==2 {next} NR>1 {input+=$2; filtered+=$3; non_chimeric+=$7} END {printf "Denoising Summary:\n - Total Input Reads: %d\n - Total Filtered Reads: %d\n - Total Non-Chimeric Reads: %d\n", input, filtered, non_chimeric}' E_rbcl_output/E_rbcl_exported_stats/stats.tsv)
loading_message "$read_summary"

# Visualize denoising stats
qiime metadata tabulate \
  --m-input-file rbcl_output/rbcl_stats.qza \
  --o-visualization rbcl_output/rbcl_stats.qzv
loading_message "Generated visualization of rbcl denoising statistics"

# Summarize feature table
qiime feature-table summarize \
  --i-table rbcl_output/rbcl_filtered_table.qza \
  --o-visualization rbcl_output/rbcl_filtered_table.qzv
loading_message "Summarized rbcl feature table"

# Export feature table
qiime tools export \
  --input-path rbcl_output/rbcl_filtered_table.qza \
  --output-path rbcl_output/rbcl_exported_table
loading_message "Exported rbcl feature table"

# Convert to TSV
biom convert \
  -i rbcl_output/rbcl_exported_table/feature-table.biom \
  -o rbcl_output/rbcl_exported_table/rbcl_feature_table.tsv \
  --to-tsv
loading_message "Converted rbcl feature table to TSV format"

# Export representative sequences
qiime tools export \
  --input-path rbcl_output/rbcl_filtered_rep_seqs.qza \
  --output-path rbcl_output/rbcl_fasta_exported
loading_message "Exported rbcl representative sequences"

# Rename FASTA
mv rbcl_output/rbcl_fasta_exported/dna-sequences.fasta rbcl_output/rbcl_fasta_exported/dada2_rbcl.fasta
loading_message "Renamed E_rbcl FASTA file to dada2_rbcl.fasta"

echo -e "\nrbcl pipeline completed successfully"
