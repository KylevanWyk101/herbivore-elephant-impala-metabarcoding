# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 10:18:33 2025
@author: kylev
"""

import pandas as pd
import matplotlib.pyplot as plt

def calculate_foo(df, level, sample_cols):
    presence = df[sample_cols].gt(0).astype(int)
    return (df[[level]].join(presence).groupby(level).max().sum(axis=1) / len(sample_cols)) * 100

def make_mirrored_plot(foo1, foo2, level, output):
    # Align and sort by max FOO descending
    taxa = sorted(set(foo1.index).union(foo2.index))
    foo1, foo2 = foo1.reindex(taxa, fill_value=0), foo2.reindex(taxa, fill_value=0)
    df = pd.DataFrame({'trnL': foo1, 'rbcL': foo2})
    df = df.assign(max=df.max(axis=1)).sort_values('max', ascending=False)
    foo1, foo2 = df['trnL'], df['rbcL']
    taxa = df.index.tolist()

    # Plot
    fig, ax = plt.subplots(figsize=(12, max(6, 0.4 * len(taxa))))  # wider plot
    y = range(len(taxa))
    ax.barh(y, -foo1, color='blue', label='trnL')
    ax.barh(y, foo2, color='red', label='rbcL')
    ax.grid(False)  # remove background grid lines

    # Add text labels
    for i, t in enumerate(taxa):
        if foo1[t] > 0:
            ax.text(-foo1[t] - 3, i, f"{foo1[t]:.0f}%", va='center', ha='right', color='blue')
        if foo2[t] > 0:
            ax.text(foo2[t] + 3, i, f"{foo2[t]:.0f}%", va='center', ha='left', color='red')

    ax.set_yticks(y)
    ax.set_yticklabels(taxa)
    ax.set_title(f'Elephant Mirrored FOO Plot by {level}')
    ax.set_xlabel('Frequency of Occurrence (%)')
    ax.axvline(0, color='black', linewidth=0.5)

    # Force symmetric axis from -100 to 100 with padding
    pad = 10
    ax.set_xlim(-120 - pad, 120 + pad)

    # Set ticks at every 25%
    ticks = range(-100, 101, 25)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(t)}%" for t in ticks])

    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.show()


# === Load data ===
trnl = pd.read_csv("E_trnl_diet_list_filtered_final.csv").fillna('NA')
rbcl = pd.read_csv("E_rbcl_diet_list_filtered_final.csv").fillna('NA')

cols_trnl = [c for c in trnl.columns if c not in ['ASV_ID', 'Species', 'Genus', 'Family', 'C3-trnl', 'Life form']]
cols_rbcl = [c for c in rbcl.columns if c not in ['ASV_ID', 'Species', 'Genus', 'Family', 'C3-rbcl', 'Life form']]

# === FOO & Plotting ===
make_mirrored_plot(
    calculate_foo(trnl, 'Family', cols_trnl),
    calculate_foo(rbcl, 'Family', cols_rbcl),
    'Family', 'E_mirrored_foo_family2.png'
)

make_mirrored_plot(
    calculate_foo(trnl, 'Genus', cols_trnl),
    calculate_foo(rbcl, 'Genus', cols_rbcl),
    'Genus', 'E_mirrored_foo_genus2.png'
)
