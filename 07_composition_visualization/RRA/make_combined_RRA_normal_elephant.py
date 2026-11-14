"""
Combined rbcL + trnL "normal" (one bar per taxon) RRA figure for Elephant --
Family level (panels A/C) and Genus/Species level (panels B/D) -- exported
directly as a print-ready PNG.

Replaces the old workflow of saving each level/gene plot separately,
screenshotting them into a 2-page PowerPoint (rbcL page, trnL page), and
exporting to PDF -- Figure_A1.pdf ended up at ~92-200 effective PPI on an
oversized 40x22.5in page despite the source PNGs being rendered at 300 DPI.
The RRA calculation logic below matches RRA_E_rbcl.py / RRA_E_trnl.py
exactly; those scripts, their CSVs, and their outputs are untouched by this
one.

Target: Wiley Ecology and Evolution appendix figure spec (<=18x22cm,
300-600 DPI, text 6-10pt). The trnL genus/species panel has 45 taxa, which
does not fit at a fully-legible 7pt/1.5x-row-pitch ratio within the height
cap -- panels with >DENSE_ROW_THRESHOLD rows drop to 6pt with tighter
spacing instead (still within Wiley's allowed range). Check the printed
mm dimensions below against 180x220mm, and inspect the dense panels'
tick labels for crowding after the first render.

Run from this directory (needs pandas + matplotlib, e.g. the
"metabarcoding" conda env on this machine):
    python make_combined_RRA_normal_elephant.py
"""
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
from PIL import Image

HERE = Path(__file__).resolve().parent

TARGET_DPI = 400
FIG_WIDTH_IN = 170 / 25.4   # under Wiley's 180mm cap
BODY_FONT_SIZE = 7          # Wiley requires figure text sized 6-10pt
HEADING_FONT_SIZE = 10
DENSE_ROW_THRESHOLD = 35

plt.rcParams.update({"font.size": BODY_FONT_SIZE})

lifeform_colors = {
    "Graminoid": "#228B22",
    "Woody taxa": "#8B4513",
    "Succulent": "#90EE90",
    "Forb": "#FF69B4",
    "NA": "#B0B0B0",
}


def get_lifeform_color(lifeform):
    return lifeform_colors.get(lifeform, "#B0B0B0")


def load_gene(csv_path):
    df = pd.read_csv(csv_path).fillna("NA")
    sample_columns = [c for c in df.columns if c not in ['ASV_ID', 'Family', 'Genus', 'Species', 'Life form']]
    total_reads = df[sample_columns].sum().sum()
    df['Genus_Species'] = df.apply(
        lambda x: f"{x['Genus']} {x['Species']}" if x['Species'] != "NA" else f"{x['Genus']}", axis=1
    )
    return df, sample_columns, total_reads


def family_rra(df, sample_columns, total_reads):
    family_lifeform = df.groupby(['Family', 'Life form'])[sample_columns].sum().reset_index()
    family_lifeform['RRA'] = (family_lifeform[sample_columns].sum(axis=1) / total_reads) * 100
    pivot_data = family_lifeform.pivot(index="Family", columns="Life form", values="RRA").fillna(0)
    pivot_data["Total"] = pivot_data.sum(axis=1)
    return pivot_data.sort_values("Total", ascending=True).drop(columns="Total")


def genus_rra(df, sample_columns, total_reads):
    genus_abundance = df.groupby(['Genus_Species', 'Life form'])[sample_columns].sum().reset_index()
    genus_abundance['RRA'] = (genus_abundance[sample_columns].sum(axis=1) / total_reads) * 100
    return genus_abundance.sort_values("RRA", ascending=True)


def panel_font_and_pitch(n_rows, dense_threshold=DENSE_ROW_THRESHOLD):
    if n_rows > dense_threshold:
        return 6, 0.105
    return 7, 0.15


def draw_family_panel(ax, pivot_data, fontsize):
    for family, row in pivot_data.iterrows():
        left = 0
        for lf in pivot_data.columns:
            value = row[lf]
            if value > 0:
                ax.barh(family, value, left=left, color=get_lifeform_color(lf), edgecolor="none")
                left += value
    ax.set_xlabel('% Relative Read Abundance', fontsize=BODY_FONT_SIZE)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.tick_params(axis='x', labelsize=BODY_FONT_SIZE)
    ax.margins(y=0.01)


def draw_genus_panel(ax, sorted_data, fontsize):
    ax.barh(
        sorted_data["Genus_Species"],
        sorted_data["RRA"],
        color=[get_lifeform_color(lf) for lf in sorted_data["Life form"]],
    )
    ax.set_xlabel('% Relative Read Abundance', fontsize=BODY_FONT_SIZE)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.tick_params(axis='x', labelsize=BODY_FONT_SIZE)
    ax.margins(y=0.01)


def main():
    rbcl_df, rbcl_cols, rbcl_total = load_gene(HERE / "E_rbcl_revised" / "E_rbcl_diet_list_filtered_final.csv")
    trnl_df, trnl_cols, trnl_total = load_gene(HERE / "E_trnl_revised" / "E_trnl_diet_list_filtered_final_fixed.csv")

    rbcl_family = family_rra(rbcl_df, rbcl_cols, rbcl_total)
    rbcl_genus = genus_rra(rbcl_df, rbcl_cols, rbcl_total)
    trnl_family = family_rra(trnl_df, trnl_cols, trnl_total)
    trnl_genus = genus_rra(trnl_df, trnl_cols, trnl_total)

    fA, pA = panel_font_and_pitch(len(rbcl_family))
    fB, pB = panel_font_and_pitch(len(rbcl_genus))
    fC, pC = panel_font_and_pitch(len(trnl_family))
    fD, pD = panel_font_and_pitch(len(trnl_genus))

    row1_h = max(pA * len(rbcl_family), pB * len(rbcl_genus), 1.2)
    row2_h = max(pC * len(trnl_family), pD * len(trnl_genus), 1.2)
    fig_height = row1_h + row2_h + 0.9  # + suptitle/legend/axis-label chrome budget

    fig = plt.figure(figsize=(FIG_WIDTH_IN, fig_height))
    gs = fig.add_gridspec(2, 2, height_ratios=[row1_h, row2_h], hspace=0.18, wspace=1.0,
                           top=0.94, bottom=0.12)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_family_panel(ax_a, rbcl_family, fA)
    draw_genus_panel(ax_b, rbcl_genus, fB)
    draw_family_panel(ax_c, trnl_family, fC)
    draw_genus_panel(ax_d, trnl_genus, fD)

    for ax, letter in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        ax.text(-0.05, 1.03, letter, transform=ax.transAxes, fontsize=HEADING_FONT_SIZE, fontweight='bold')

    ax_a.set_title(r"$\mathit{rbc}$L — Family", fontsize=BODY_FONT_SIZE)
    ax_b.set_title(r"$\mathit{rbc}$L — Genus/Species", fontsize=BODY_FONT_SIZE)
    ax_c.set_title(r"$\mathit{trn}$L — Family", fontsize=BODY_FONT_SIZE)
    ax_d.set_title(r"$\mathit{trn}$L — Genus/Species", fontsize=BODY_FONT_SIZE)

    fig.suptitle("Figure A1", fontsize=HEADING_FONT_SIZE, y=0.99)

    present = (set(rbcl_family.columns) | set(trnl_family.columns) |
               set(rbcl_genus["Life form"]) | set(trnl_genus["Life form"]))
    present_lifeforms = [lf for lf in lifeform_colors if lf in present]
    handles = [mpatches.Patch(color=lifeform_colors[lf], label=lf) for lf in present_lifeforms]
    fig.legend(handles=handles, title="Life form", loc="lower center",
               ncol=len(handles), bbox_to_anchor=(0.5, 0.01), fontsize=BODY_FONT_SIZE,
               title_fontsize=BODY_FONT_SIZE, frameon=False)

    out_path = HERE / "Figure_A1_composite.png"
    fig.savefig(out_path, dpi=TARGET_DPI, bbox_inches='tight')

    with Image.open(out_path) as im:
        w_mm, h_mm = im.width / TARGET_DPI * 25.4, im.height / TARGET_DPI * 25.4
        size_kb = out_path.stat().st_size / 1024
        print(f"Saved {out_path}")
        print(f"  {im.width}x{im.height}px @ {TARGET_DPI} DPI -> {w_mm:.1f} x {h_mm:.1f} mm, {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
