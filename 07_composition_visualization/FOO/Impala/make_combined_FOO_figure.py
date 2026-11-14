"""
Combined Family (panel A) + Genus (panel B) mirrored FOO bar chart for
Impala, exported directly as a print-ready PNG.

Replaces the old workflow of saving each level's plot separately, screenshotting
them into PowerPoint side by side, and exporting to PDF -- Figure_A6.pdf ended up
at ~95-165 effective PPI despite the source PNGs being rendered at 600 DPI, the
same PowerPoint-export resolution loss already found and fixed for the Sankey
figures in Statistics/niche_partitioning/. The FOO calculation and mirrored-bar
logic below matches FOO_calc_plot.py exactly; that script and its outputs are
untouched by this one.

Run from this directory (needs pandas + matplotlib + Pillow, e.g. the "humann"
conda env on this machine):
    python make_combined_FOO_figure.py
"""
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

TARGET_DPI = 400              # Wiley Ecology and Evolution: 300-600 DPI
FIG_WIDTH_IN = 170 / 25.4     # Wiley appendix figure cap is 18cm wide; stay under it
MAX_HEIGHT_IN = 200 / 25.4    # cap is 22cm tall; leave a margin for bbox_inches='tight' growth
ROW_HEIGHT_IN = 0.15          # vertical space per taxon row, sized for BODY_FONT_SIZE below
BODY_FONT_SIZE = 7            # Wiley requires figure text sized 6-10pt at final print size
HEADING_FONT_SIZE = 10        # figure title / panel letters, top of Wiley's allowed range

plt.rcParams.update({"font.size": BODY_FONT_SIZE})


def calculate_foo(df, level, sample_cols):
    presence = df[sample_cols].gt(0).astype(int)
    return (df[[level]].join(presence).groupby(level).max().sum(axis=1) / len(sample_cols)) * 100


def draw_panel(ax, foo1, foo2):
    """Draw one mirrored trnL/rbcL FOO bar chart onto an existing Axes.
    Same logic as make_mirrored_plot() in FOO_calc_plot.py, minus the
    per-panel title (dropped in favour of the shared "A"/"B" labels) and
    the figure creation/save (this is now one panel of a shared figure)."""
    taxa = sorted(set(foo1.index).union(foo2.index))
    foo1, foo2 = foo1.reindex(taxa, fill_value=0), foo2.reindex(taxa, fill_value=0)
    df = pd.DataFrame({'trnL': foo1, 'rbcL': foo2})
    df = df.assign(max=df.max(axis=1)).sort_values('max', ascending=False)
    foo1, foo2 = df['trnL'], df['rbcL']
    taxa = df.index.tolist()

    y = range(len(taxa))
    ax.barh(y, -foo1, color='blue', label=r'$\mathit{trn}$L')
    ax.barh(y, foo2, color='red', label=r'$\mathit{rbc}$L')
    ax.grid(False)

    for i, t in enumerate(taxa):
        if foo1[t] > 0:
            ax.text(-foo1[t] - 3, i, f"{foo1[t]:.0f}%", va='center', ha='right', color='blue')
        if foo2[t] > 0:
            ax.text(foo2[t] + 3, i, f"{foo2[t]:.0f}%", va='center', ha='left', color='red')

    ax.set_yticks(y)
    ax.set_yticklabels(taxa)
    ax.set_xlabel('Frequency of Occurrence (%)')
    ax.axvline(0, color='black', linewidth=0.5)

    pad = 20
    ax.set_xlim(-120 - pad, 120 + pad)
    ticks = range(-100, 101, 50)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(t)}%" for t in ticks])
    ax.legend(loc='upper right')

    return len(taxa)


def main():
    trnl = pd.read_csv(HERE / "I_trnl_diet_list_filtered_final.csv").fillna('NA')
    rbcl = pd.read_csv(HERE / "I_rbcl_diet_list_filtered_final.csv").fillna('NA')

    cols_trnl = [c for c in trnl.columns if c not in ['ASV_ID', 'Species', 'Genus', 'Family', 'C1-trnl', '48-trnl', 'Life form']]
    cols_rbcl = [c for c in rbcl.columns if c not in ['ASV_ID', 'Species', 'Genus', 'Family', 'C1-rbcl', 'Life form']]

    foo1_family = calculate_foo(trnl, 'Family', cols_trnl)
    foo2_family = calculate_foo(rbcl, 'Family', cols_rbcl)
    foo1_genus = calculate_foo(trnl, 'Genus', cols_trnl)
    foo2_genus = calculate_foo(rbcl, 'Genus', cols_rbcl)

    n_genus_taxa = len(set(foo1_genus.index) | set(foo2_genus.index))
    fig_height = min(MAX_HEIGHT_IN, max(3, ROW_HEIGHT_IN * n_genus_taxa))  # sized for the denser genus-level panel, capped to fit the page

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(FIG_WIDTH_IN, fig_height))

    draw_panel(ax_a, foo1_family, foo2_family)
    draw_panel(ax_b, foo1_genus, foo2_genus)

    fig.suptitle("Figure A6", fontsize=HEADING_FONT_SIZE, y=1.0)
    ax_a.text(-0.05, 1.02, "A", transform=ax_a.transAxes, fontsize=HEADING_FONT_SIZE, fontweight='bold')
    ax_b.text(-0.05, 1.02, "B", transform=ax_b.transAxes, fontsize=HEADING_FONT_SIZE, fontweight='bold')

    plt.tight_layout()

    out_path = HERE / "Figure_A6_composite.png"
    fig.savefig(out_path, dpi=TARGET_DPI, bbox_inches='tight')

    from PIL import Image
    with Image.open(out_path) as im:
        w_in, h_in = im.width / TARGET_DPI, im.height / TARGET_DPI
        size_kb = out_path.stat().st_size / 1024
        print(f"Saved {out_path}")
        print(f"  {im.width}x{im.height}px @ {TARGET_DPI} DPI -> "
              f"{w_in*25.4:.1f} x {h_in*25.4:.1f} mm, {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
