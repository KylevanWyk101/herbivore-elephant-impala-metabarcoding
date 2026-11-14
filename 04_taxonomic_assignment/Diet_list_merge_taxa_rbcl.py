import pandas as pd

# Load files
hits = pd.read_csv("rbcl_final_hits.csv")
feature = pd.read_csv("asv_table_filtered_rbcl.tsv", sep="\t")
feature.columns = ["ASV_ID"] + list(feature.columns[1:])  # Rename first column to ASV_ID

# Merge taxonomy with feature table
df = feature.merge(hits, on="ASV_ID", how="left")

# Fill missing values in hits columns with "NA"
hits_columns = [col for col in hits.columns if col != "ASV_ID"]
df[hits_columns] = df[hits_columns].fillna("NA")

# Define sample sets
elephant_samples = [str(i) + "-rbcl" for i in range(49, 73)] + ['C3-rbcl']
impala_samples = [str(i) + "-rbcl" for i in range(25, 49)] + ['C1-rbcl']

# Filter sample lists to include only columns that exist in the DataFrame
elephant_samples = [col for col in elephant_samples if col in df.columns]
impala_samples = [col for col in impala_samples if col in df.columns]

# Filter ASVs with non-zero abundance for each animal group
elephant_df = df[['ASV_ID'] + hits_columns + elephant_samples]
elephant_df = elephant_df[elephant_df[elephant_samples].sum(axis=1) > 0]

impala_df = df[['ASV_ID'] + hits_columns + impala_samples]
impala_df = impala_df[impala_df[impala_samples].sum(axis=1) > 0]

# Save outputs
df.to_csv("full_rbcl_diet_list.csv", index=False)
elephant_df.to_csv("E_rbcl_diet_list.csv", index=False)
impala_df.to_csv("I_rbcl_diet_list.csv", index=False)

print("✅ Files saved: full_rbcl_diet_list.csv, E_rbcl_diet_list.csv, I_rbcl_diet_list.csv")
