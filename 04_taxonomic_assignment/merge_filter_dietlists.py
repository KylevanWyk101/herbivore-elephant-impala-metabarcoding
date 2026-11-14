import pandas as pd
###  RBCL  ###
# Load the CSV files
impala_df = pd.read_csv('I_rbcl_diet_list_filtered_final.csv')
elephant_df = pd.read_csv('E_rbcl_diet_list_filtered_final.csv')

# Define taxonomy columns
taxonomy_cols = ['ASV_ID', 'Family', 'Genus', 'Species', 'Life form']

# Remove rows where Genus is blank or NA in each dataset and report removed ASVs
impala_df_original = impala_df.copy()
impala_df = impala_df[impala_df['Genus'].notna() & (impala_df['Genus'] != '')]
removed_impala_asvs = impala_df_original[~impala_df_original['ASV_ID'].isin(impala_df['ASV_ID'])]['ASV_ID'].tolist()
print(f"Removed {len(removed_impala_asvs)} rows from impala dataset due to blank or NA Genus:")
if removed_impala_asvs:
    print(f"Removed ASVs (impala): {', '.join(map(str, removed_impala_asvs))}")
else:
    print("No ASVs removed from impala dataset.")

elephant_df_original = elephant_df.copy()
elephant_df = elephant_df[elephant_df['Genus'].notna() & (elephant_df['Genus'] != '')]
removed_elephant_asvs = elephant_df_original[~elephant_df_original['ASV_ID'].isin(elephant_df['ASV_ID'])]['ASV_ID'].tolist()
print(f"Removed {len(removed_elephant_asvs)} rows from elephant dataset due to blank or NA Genus:")
if removed_elephant_asvs:
    print(f"Removed ASVs (elephant): {', '.join(map(str, removed_elephant_asvs))}")
else:
    print("No ASVs removed from elephant dataset.")

# Identify sample columns (excluding taxonomy columns)
impala_sample_cols = [col for col in impala_df.columns if col not in taxonomy_cols]
elephant_sample_cols = [col for col in elephant_df.columns if col not in taxonomy_cols]

# Rename sample columns to include _I for impala and _E for elephant
impala_df = impala_df.rename(columns={col: f"{col}_I" for col in impala_sample_cols})
elephant_df = elephant_df.rename(columns={col: f"{col}_E" for col in elephant_sample_cols})

# Merge the datasets on taxonomy columns using an outer join
merged_df = pd.merge(
    impala_df,
    elephant_df,
    on=taxonomy_cols,
    how='outer'
)

# Fill missing values in sample columns with 0
merged_df = merged_df.fillna(0)

# Remove columns (excluding taxonomy columns) that contain only zeros and report them
sample_cols = [col for col in merged_df.columns if col not in taxonomy_cols]
zero_cols = [col for col in sample_cols if merged_df[col].sum() == 0]
non_zero_cols = [col for col in sample_cols if col not in zero_cols]
# Keep taxonomy columns and non-zero sample columns
final_cols = taxonomy_cols + non_zero_cols
merged_df = merged_df[final_cols]

print(f"Removed {len(zero_cols)} sample columns with only zeros:")
if zero_cols:
    print(f"Removed columns: {', '.join(zero_cols)}")
else:
    print("No sample columns removed.")

# Save the combined dataset to a new CSV
merged_df.to_csv('rbcl_combined_dietlist.csv', index=False)
print("Combined dataset saved as 'rbcl.csv'")



#### TRNL #####
# Load the CSV files
impala_df = pd.read_csv('I_trnl_diet_list_filtered_final.csv')
elephant_df = pd.read_csv('E_trnl_diet_list_filtered_final.csv')

# Define taxonomy columns
taxonomy_cols = ['ASV_ID', 'Family', 'Genus', 'Species', 'Life form']

# Remove rows where Genus is blank or NA in each dataset and report removed ASVs
impala_df_original = impala_df.copy()
impala_df = impala_df[impala_df['Genus'].notna() & (impala_df['Genus'] != '')]
removed_impala_asvs = impala_df_original[~impala_df_original['ASV_ID'].isin(impala_df['ASV_ID'])]['ASV_ID'].tolist()
print(f"Removed {len(removed_impala_asvs)} rows from impala dataset due to blank or NA Genus:")
if removed_impala_asvs:
    print(f"Removed ASVs (impala): {', '.join(map(str, removed_impala_asvs))}")
else:
    print("No ASVs removed from impala dataset.")

elephant_df_original = elephant_df.copy()
elephant_df = elephant_df[elephant_df['Genus'].notna() & (elephant_df['Genus'] != '')]
removed_elephant_asvs = elephant_df_original[~elephant_df_original['ASV_ID'].isin(elephant_df['ASV_ID'])]['ASV_ID'].tolist()
print(f"Removed {len(removed_elephant_asvs)} rows from elephant dataset due to blank or NA Genus:")
if removed_elephant_asvs:
    print(f"Removed ASVs (elephant): {', '.join(map(str, removed_elephant_asvs))}")
else:
    print("No ASVs removed from elephant dataset.")

# Identify sample columns (excluding taxonomy columns)
impala_sample_cols = [col for col in impala_df.columns if col not in taxonomy_cols]
elephant_sample_cols = [col for col in elephant_df.columns if col not in taxonomy_cols]

# Rename sample columns to include _I for impala and _E for elephant
impala_df = impala_df.rename(columns={col: f"{col}_I" for col in impala_sample_cols})
elephant_df = elephant_df.rename(columns={col: f"{col}_E" for col in elephant_sample_cols})

# Merge the datasets on taxonomy columns using an outer join
merged_df = pd.merge(
    impala_df,
    elephant_df,
    on=taxonomy_cols,
    how='outer'
)

# Fill missing values in sample columns with 0
merged_df = merged_df.fillna(0)

# Remove columns (excluding taxonomy columns) that contain only zeros and report them
sample_cols = [col for col in merged_df.columns if col not in taxonomy_cols]
zero_cols = [col for col in sample_cols if merged_df[col].sum() == 0]
non_zero_cols = [col for col in sample_cols if col not in zero_cols]
# Keep taxonomy columns and non-zero sample columns
final_cols = taxonomy_cols + non_zero_cols
merged_df = merged_df[final_cols]

print(f"Removed {len(zero_cols)} sample columns with only zeros:")
if zero_cols:
    print(f"Removed columns: {', '.join(zero_cols)}")
else:
    print("No sample columns removed.")

# Save the combined dataset to a new CSV
merged_df.to_csv('trnl_combined_dietlist.csv', index=False)
print("Combined dataset saved as 'trnl_combined_dietlist.csv'")