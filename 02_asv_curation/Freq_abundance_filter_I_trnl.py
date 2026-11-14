# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 15:33:20 2025

@author: kylev
"""

import pandas as pd

# Load ASV table
asv_table = pd.read_csv("I_trnl_LULU_curated_table.tsv", sep="\t", index_col=0)

# -------- Step 1: Frequency Filter --------
# Remove ASVs that appear in only one sample (i.e., only one non-zero value per row)
asv_filtered = asv_table[asv_table.astype(bool).sum(axis=1) > 1]


# -------- Save the filtered ASV table --------
asv_filtered.to_csv("asv_table_filtered_I_trnl.tsv", sep="\t")