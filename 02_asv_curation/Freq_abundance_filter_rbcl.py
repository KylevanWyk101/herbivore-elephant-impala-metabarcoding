# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 15:33:20 2025

@author: kylev
"""

import pandas as pd

# Load ASV table
asv_table = pd.read_csv("rbcl_LULU_curated_table.tsv", sep="\t", index_col=0)

# -------- Step 1: Frequency Filter --------
# Remove ASVs that appear in only one sample (i.e., only one non-zero value per row)
asv_filtered = asv_table[asv_table.astype(bool).sum(axis=1) > 1]


# -------- Step 2: Abundance Filter (rbcL only) --------
# Calculate total reads across all samples
total_reads = asv_filtered.sum().sum()

# Filter threshold: 0.005% of total reads
threshold = 0.00005 * total_reads

# Remove ASVs whose total read count across all samples is below the threshold
asv_final = asv_filtered[asv_filtered.sum(axis=1) >= threshold]

# -------- Save the filtered ASV table --------
asv_final.to_csv("asv_table_filtered_rbcl.tsv", sep="\t")

asv_table_final = pd.read_csv("asv_table_filtered_rbcl.tsv", sep="\t", index_col=0)
