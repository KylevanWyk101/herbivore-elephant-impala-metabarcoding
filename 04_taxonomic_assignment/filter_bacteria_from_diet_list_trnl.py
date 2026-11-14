import pandas as pd
# Now filter the diet list to remove bacterial assignments

# ELEPHANT DATA
# File paths
taxonomy_file = "trnl_blastn_genbank_filtered_results.csv"
diet_file = "E_trnl_diet_list.csv"
output_file = "E_trnl_diet_list_filtered.csv"

# Read taxonomy file and get unique ASV IDs
taxonomy_df = pd.read_csv(taxonomy_file)
bacterial_asvs = taxonomy_df['query id'].unique().tolist()
print(f"Unique ASVs in taxonomy file: {len(bacterial_asvs)}")

# Read dietary TSV and filter out bacterial ASVs
diet_df = pd.read_csv(diet_file)
filtered_diet_df = diet_df[~diet_df['ASV_ID'].isin(bacterial_asvs)]

# Report counts
removed_count = len(diet_df) - len(filtered_diet_df)
print(f"Removed {removed_count} bacterial ASVs from dietary CSV")
print(f"Original rows: {len(diet_df)}, Filtered rows: {len(filtered_diet_df)}")

# Save filtered CSV
filtered_diet_df.to_csv(output_file, index=False)
print(f"Filtered dietary CSV saved to {output_file}")



# IMPALA data

# File paths
taxonomy_file = "trnl_blastn_genbank_filtered_results.csv"
diet_file = "I_trnl_diet_list.csv"
output_file = "I_trnl_diet_list_filtered.csv"

# Read taxonomy file and get unique ASV IDs
taxonomy_df = pd.read_csv(taxonomy_file)
bacterial_asvs = taxonomy_df['query id'].unique().tolist()
print(f"Unique ASVs in taxonomy file: {len(bacterial_asvs)}")

# Read dietary TSV and filter out bacterial ASVs
diet_df = pd.read_csv(diet_file)
filtered_diet_df = diet_df[~diet_df['ASV_ID'].isin(bacterial_asvs)]

# Report counts
removed_count = len(diet_df) - len(filtered_diet_df)
print(f"Removed {removed_count} bacterial ASVs from dietary CSV")
print(f"Original rows: {len(diet_df)}, Filtered rows: {len(filtered_diet_df)}")

# Save filtered CSV
filtered_diet_df.to_csv(output_file, index=False)
print(f"Filtered dietary CSV saved to {output_file}")