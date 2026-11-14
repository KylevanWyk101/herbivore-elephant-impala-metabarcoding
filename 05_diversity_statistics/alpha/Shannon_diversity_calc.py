# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 14:47:51 2025

@author: kylev
"""
import pandas as pd
from skbio.diversity import alpha_diversity

# List of RRA tables and their metadata
datasets = [
    {'file': 'E_rbcl_diet_rra_only.csv', 'group': 'Elephant', 'locus': 'rbcL'},
    {'file': 'I_rbcl_diet_rra_only.csv', 'group': 'Impala', 'locus': 'rbcL'},
    {'file': 'E_trnl_diet_rra_only.csv', 'group': 'Elephant', 'locus': 'trnL'},
    {'file': 'I_trnl_diet_rra_only.csv', 'group': 'Impala', 'locus': 'trnL'}
]

# Initialize list to store Shannon index results
all_shannon_results = []

# Process each dataset
for dataset in datasets:
    # Load the RRA table
    df = pd.read_csv(dataset['file'])
    
    # Identify RRA columns (those ending with '_RRA')
    rra_columns = [col for col in df.columns if col.endswith('_RRA')]
    
  # Calculate Shannon index at ASV level
    shannon_index = alpha_diversity('shannon', df[rra_columns].values.T, ids=rra_columns)

    # Create results dataframe
    sample_ids = [col.replace('_RRA', '') for col in rra_columns]
    shannon_results = pd.DataFrame({
        'SampleID': sample_ids,
        'ShannonIndex': shannon_index,
        'Group': dataset['group'],
        'Locus': dataset['locus']
    })
    
    # Append to list
    all_shannon_results.append(shannon_results)

# Combine all Shannon index results
combined_shannon = pd.concat(all_shannon_results, ignore_index=True)

# Save combined results to a CSV file
combined_shannon.to_csv('combined_shannon_index.csv', index=False)

print("Shannon index calculation complete for all datasets.")
print("Combined results saved to 'combined_shannon_index_.csv'.")
print(combined_shannon)
