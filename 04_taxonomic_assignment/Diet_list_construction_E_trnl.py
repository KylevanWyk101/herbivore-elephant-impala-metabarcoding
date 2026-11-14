import pandas as pd

# Load the feature table (ASV IDs in first column)
feature_table = pd.read_csv("asv_table_filtered_E_trnl.tsv", sep="\t", index_col=0)

# Load taxonomy assignments and rename ID column for merging
taxonomy_df = pd.read_csv("asv_taxonomy_assignments.csv").rename(columns={"Feature ID": "ASV_ID"})

# Merge on ASV ID (ASV_ID from taxonomy_df, index from feature_table)
merged_df = feature_table.merge(taxonomy_df, how="left", left_index=True, right_on="ASV_ID")

# Fill NaNs in taxonomy columns with 'unclassified'
taxonomy_cols = ['Species', 'Genus', 'Family', 'Life form']
merged_df[taxonomy_cols] = merged_df[taxonomy_cols].fillna('unclassified')

# Reorder columns: taxonomy first, then sample columns
sample_columns = [col for col in feature_table.columns]
ordered_columns = ['ASV_ID', 'Species', 'Genus', 'Family', 'Life form'] + sample_columns
merged_df = merged_df[ordered_columns]

# Save final table
merged_df.to_csv("E_trnl_diet_list.csv", index=False)