import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
import numpy as np

# =====================================================
# === LOAD DATA ===
# =====================================================
df = pd.read_csv("I_rbcl_diet_list_filtered_final.csv").fillna("NA")

# Identify sample columns
sample_columns = [col for col in df.columns if col not in ['ASV_ID','Family','Genus','Species','Life form']]

# Total reads across all samples
total_reads = df[sample_columns].sum().sum()

# =====================================================
# === LIFE FORM COLOR MAP ===
# =====================================================
lifeform_colors = {
    "Graminoid": "#228B22",   # dark green
    "Woody taxa": "#8B4513",  # brown
    "Succulent": "#90EE90",   # light green
    "Forb": "#FF69B4",        # pink
    "NA": "#B0B0B0"           # gray
}

def get_lifeform_color(lifeform):
    return lifeform_colors.get(lifeform, "#B0B0B0")

def build_lifeform_legend(present_lifeforms):
    """Build legend handles only for life forms present in the dataset."""
    return [mpatches.Patch(color=lifeform_colors[lf], label=lf) for lf in present_lifeforms]

# =====================================================
# === RELATIVE READ ABUNDANCE BY FAMILY (TIGHT BARS) ===
# =====================================================
family_lifeform = df.groupby(['Family', 'Life form'])[sample_columns].sum().reset_index()
family_lifeform['RRA'] = (family_lifeform[sample_columns].sum(axis=1) / total_reads) * 100

# --- Save to CSV ---
family_rra_table = family_lifeform[['Family', 'Life form', 'RRA']].sort_values('RRA', ascending=False)
family_rra_table.to_csv("Impala_rbcl_RRA_by_family.csv", index=False)

# --- Plot ---
pivot_data = family_lifeform.pivot(index="Family", columns="Life form", values="RRA").fillna(0)
pivot_data["Total"] = pivot_data.sum(axis=1)
pivot_data = pivot_data.sort_values("Total", ascending=True).drop(columns="Total")

n_families = pivot_data.shape[0]
fig_height = max(5, n_families * 0.3)

plt.figure(figsize=(9, fig_height))
for i, (family, row) in enumerate(pivot_data.iterrows()):
    left = 0
    for lf in pivot_data.columns:
        value = row[lf]
        if value > 0:
            plt.barh(family, value, left=left, color=get_lifeform_color(lf), edgecolor="none")
            left += value

plt.title(r'Impala $\mathit{rbcl}$ Total RRA by Family')
plt.xlabel('% Relative Read Abundance')
plt.ylabel('Family')
plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
present_lifeforms = [lf for lf in pivot_data.columns if lf in lifeform_colors]
plt.legend(handles=build_lifeform_legend(present_lifeforms),
           title="Life form", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.margins(y=0.01)
# plt.savefig("I_rbcl_rra_by_family_stacked.png", dpi=300)
plt.show()

# =====================================================
# === RELATIVE READ ABUNDANCE BY GENUS (WITH SPECIES) ===
# =====================================================
df['Genus_Species'] = df.apply(
    lambda x: f"{x['Genus']} {x['Species']}" if x['Species'] != "NA" else f"{x['Genus']}", axis=1
)

genus_abundance = df.groupby(['Genus_Species','Life form'])[sample_columns].sum().reset_index()
genus_abundance['RRA'] = (genus_abundance[sample_columns].sum(axis=1) / total_reads) * 100

# --- Save to CSV ---
genus_rra_table = genus_abundance[['Genus_Species', 'Life form', 'RRA']].sort_values('RRA', ascending=False)
genus_rra_table.to_csv("Impala_rbcl_RRA_by_genus_species.csv", index=False)

# --- Plot ---
sorted_data = genus_abundance.sort_values("RRA", ascending=True)
n_taxa = sorted_data["Genus_Species"].nunique()
fig_height = max(5, n_taxa * 0.3)

plt.figure(figsize=(9, fig_height))
plt.barh(
    sorted_data["Genus_Species"],
    sorted_data["RRA"],
    color=[get_lifeform_color(lf) for lf in sorted_data["Life form"]]
)

plt.title(r'Impala $\mathit{rbcl}$ Total RRA by Genus Species')
plt.xlabel('% Relative Read Abundance')
plt.ylabel('Genus Species')
plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
present_lifeforms = sorted_data["Life form"].unique()
plt.legend(handles=build_lifeform_legend(present_lifeforms),
           title="Life form", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.margins(y=0.01)
# plt.savefig("I_rbcl_rra_by_genus_species.png", dpi=300)
plt.show()

# =====================================================
# === COLOR PALETTE (tab20 + tab20b + tab20c, no greys) ===
# =====================================================
cmap_20  = plt.cm.tab20.colors
cmap_b   = plt.cm.tab20b.colors
cmap_c   = plt.cm.tab20c.colors
all_colors = list(cmap_20) + list(cmap_b) + list(cmap_c)

def is_gray(rgb, tol=0.1):
    return max(rgb) - min(rgb) < tol

filtered_colors = [c for c in all_colors if not is_gray(c)]

# =====================================================
# === STACKED BARPLOT: FAMILY RRA PER SAMPLE ===
# =====================================================
family_sample = df.groupby('Family')[sample_columns].sum().T
family_sample_percent = family_sample.div(family_sample.sum(axis=1), axis=0) * 100
families = family_sample_percent.columns

color_map_family = {fam: ("#B0B0B0" if fam == "NA" else filtered_colors[i % len(filtered_colors)])
                    for i, fam in enumerate(families)}

fig, ax = plt.subplots(figsize=(16, 9))
family_sample_percent.plot(
    kind='bar', stacked=True, ax=ax,
    color=[color_map_family[f] for f in families]
)

plt.title(r'Impala $\mathit{rbcl}$ RRA per Sample by Family', fontsize=14)
ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('% Relative Read Abundance', fontsize=12)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
plt.xticks(rotation=45, ha='right')

plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, -0.15),
           fontsize=12, title_fontsize=10, ncol=6, frameon=False)
plt.tight_layout()
plt.savefig("stacked_rra_per_sample_by_family.png", dpi=300)
plt.show()

# =====================================================
# === STACKED BARPLOT: GENUS (SPECIES) RRA PER SAMPLE ===
# =====================================================
genus_sample = df.groupby('Genus_Species')[sample_columns].sum().T
genus_sample_percent = genus_sample.div(genus_sample.sum(axis=1), axis=0) * 100

genera = [g.split()[0] if g != "NA" else "NA" for g in genus_sample_percent.columns]
genus_species_labels = genus_sample_percent.columns

unique_genera = sorted(set([g for g in genera if g != "NA"]))
color_map_genus_species = {}

for i, genus in enumerate(unique_genera):
    base_color = filtered_colors[i % len(filtered_colors)]
    for label in genus_species_labels:
        if label.startswith(genus + " "):
            color_map_genus_species[label] = base_color

for label in genus_species_labels:
    if label == "NA":
        color_map_genus_species[label] = "#B0B0B0"
    elif label not in color_map_genus_species:
        color_map_genus_species[label] = filtered_colors[hash(label) % len(filtered_colors)]

fig, ax = plt.subplots(figsize=(18, 10))
genus_sample_percent.plot(
    kind='bar', stacked=True, ax=ax,
    color=[color_map_genus_species[g] for g in genus_species_labels]
)

plt.title(r'Impala $\mathit{rbcl}$ RRA per Sample by Genus Species', fontsize=14)
ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('% Relative Read Abundance', fontsize=12)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
plt.xticks(rotation=45, ha='right')

plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, -0.15),
           fontsize=12, title_fontsize=10, ncol=6, frameon=False)

plt.tight_layout()
plt.savefig("stacked_rra_per_sample_by_genus_species.png", dpi=300)
plt.show()
