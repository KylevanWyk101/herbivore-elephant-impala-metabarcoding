# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 13:31:47 2025

@author: kylev
"""
from Bio import SeqIO

input_fasta = "dada2_I_trnl.fasta"
output_fasta = "dada2_I_trnl_min50.fasta"
min_len = 50

total = 0
kept = 0
discarded = 0

with open(output_fasta, "w") as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        total += 1
        if len(record.seq) >= min_len:
            kept += 1
            SeqIO.write(record, out_handle, "fasta")
        else:
            discarded += 1

print("📊 Length-based filtering summary:")
print(f" - Total input sequences:      {total}")
print(f" - Sequences retained (≥{min_len} bp): {kept}")
print(f" - Sequences discarded (<{min_len} bp): {discarded}")