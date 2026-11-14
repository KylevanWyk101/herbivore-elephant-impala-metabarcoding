import pandas as pd

# Load the TSV file
df = pd.read_csv('E_trnl_diet_list_filtered_final.csv').fillna('NA')

# Step 1: Remove rows where Genus is "unclassified" (contaminants)
df = df[~df['Genus'].eq('NA')]

# Step 2: Identify sample columns (those after 'Family' but before the last column, assumed control)
family_index = df.columns.get_loc('Life form')
sample_columns = df.columns[family_index + 1:-1].tolist()

# Step 3: Recalculate total reads per sample after filtering
total_reads = df[sample_columns].sum()

# Step 4: Calculate RRA (Relative Read Abundance) for each sample
for col in sample_columns:
    df[f"{col}_RRA"] = df[col] / total_reads[col]

# Step 5: Drop RRA columns that are entirely NA
all_rra_columns = [f"{col}_RRA" for col in sample_columns]
rra_df = df[all_rra_columns]
na_rra_columns = rra_df.columns[rra_df.isna().all()].tolist()
df.drop(columns=na_rra_columns, inplace=True)

# Step 6: Retain only taxonomy + RRA columns
taxonomic_columns = ['ASV_ID', 'Species', 'Genus', 'Family','Life form']
retained_rra_columns = [col for col in all_rra_columns if col not in na_rra_columns]
df = df[taxonomic_columns + retained_rra_columns]

# Step 7: Save final output
df.to_csv('E_trnl_diet_rra_only.csv', index=False)

# Step 8: Report dropped samples
dropped_sample_names = [col.replace('_RRA', '') for col in na_rra_columns]

# Final message
print("Processing complete. Output saved to 'E_trnl_diet_rra_only.csv'.")
print(f"Processed {len(retained_rra_columns)} RRA columns.")
print(f"Dropped {len(dropped_sample_names)} samples due to all NA values: {', '.join(dropped_sample_names)}")

