import pandas as pd

# Load the ASV summary table
summary_df = pd.read_csv('ASV_summary_table_rbcL.csv')

# Filter rows where the category contains '(gain)'
gain_df = summary_df[summary_df['Category'].str.contains('\(gain\)', na=False)]

# Load the elephant diet list
elephant_df = pd.read_csv('elephant_rbcl_diet_list.csv')
elephant_asvs = set(elephant_df['ASV_ID'])

# Load the impala diet list
impala_df = pd.read_csv('impala_rbcl_diet_list.csv')
impala_asvs = set(impala_df['ASV_ID'])

# Create a list to hold the results
results = []

# For each gain ASV, check presence in diet lists and extract rows from elephant diet list
for asv in gain_df['ASV_ID']:
    # Extract all rows from elephant diet list for this ASV
    ele_rows = elephant_df[elephant_df['ASV_ID'] == asv].copy()
    if asv in elephant_asvs and asv in impala_asvs:
        category = 'Both'
    elif asv in elephant_asvs:
        category = 'Elephant'
    elif asv in impala_asvs:
        category = 'Impala'
    else:
        category = 'Neither'
    if not ele_rows.empty:
        ele_rows['Diet Presence'] = category
        results.extend(ele_rows.to_dict('records'))
    else:
        results.append({'ASV_ID': asv, 'Diet Presence': category})

# Create a new DataFrame from results
new_df = pd.DataFrame(results)

# Save to new CSV
new_df.to_csv('gain_asvs_diet_presence_elephant_details.csv', index=False)

print("New CSV file 'gain_asvs_diet_presence_elephant_details.csv' created successfully!")