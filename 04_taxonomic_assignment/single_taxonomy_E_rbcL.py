import pandas as pd

# Load the data
df = pd.read_csv("E_rbcl_final_hits.csv")

# Check that required columns exist
required_columns = ['query id', 'Species', 'Genus', 'Family', 'Life form', 'confidence score']
if not all(col in df.columns for col in required_columns):
    raise ValueError("Missing required columns in the input file.")

# Get max confidence score rows
max_conf = df.groupby('query id')['confidence score'].max().reset_index()
df_max = df.merge(max_conf, on=['query id', 'confidence score'])

# Function to resolve taxonomy based on top hits
def resolve_taxonomy(group, query_id):
    result = {'Feature ID': query_id}

    if group['Family'].nunique() == 1:
        result['Family'] = group['Family'].iloc[0]
        if group['Genus'].nunique() == 1:
            result['Genus'] = group['Genus'].iloc[0]
            if group['Species'].nunique() == 1:
                result['Species'] = group['Species'].iloc[0]
            else:
                result['Species'] = 'unknown'
        else:
            result['Genus'] = 'unknown'
            result['Species'] = 'unknown'
        # Life form is safe to copy if Family is consistent
        result['Life form'] = group['Life form'].iloc[0]
    else:
        # Conflict already at Family level
        result['Family'] = 'unknown'
        result['Genus'] = 'unknown'
        result['Species'] = 'unknown'
        result['Life form'] = 'unknown'

    return pd.Series(result)

# Preserve original ASV order
original_order = df['query id'].drop_duplicates()

# Apply taxonomy resolution
taxonomy_df = (
    df_max.groupby('query id', group_keys=False)
    .apply(lambda g: resolve_taxonomy(g.drop(columns='query id'), g.name))
    .reset_index(drop=True)
)

# Restore original order
taxonomy_df['Feature ID'] = pd.Categorical(taxonomy_df['Feature ID'], categories=original_order, ordered=True)
taxonomy_df = taxonomy_df.sort_values('Feature ID')

# Save output
taxonomy_df.to_csv("asv_taxonomy_assignments2.csv", index=False)
