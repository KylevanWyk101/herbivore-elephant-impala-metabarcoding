"""
Combined rbcL + trnL per-sample stacked RRA figure for Impala -- Family
level (panels A/C) and Genus/Species level (panels B/D) -- exported
directly as a print-ready PNG.

Replaces the old workflow of saving each level/gene plot separately,
screenshotting them into a 2-page PowerPoint (rbcL page, trnL page), and
exporting to PDF -- Figure_A4.pdf ended up at ~203-238 effective PPI on an
oversized 40x22.5in page despite the source PNGs being rendered at 300 DPI.
The per-sample normalization and colour-palette logic below matches
RRA_I_rbcl.py / RRA_I_trnl.py exactly (every taxon keeps its own
distinguishable colour from the same tab20+tab20b+tab20c palette, no data
or grouping changes); those scripts, their CSVs, and their outputs are
untouched by this one.

A shared gradient-by-rank colour scheme was tried and rejected -- adjacent
ranks looked too similar to tell segments apart, defeating the point of a
stacked bar. Back to full per-taxon distinguishable colours + full legend,
accepting that this figure needs most/all of the 22cm height budget given
it has to legend up to ~43 individually-coloured taxa.

Deliberate deviations from the originals (see project plan for why):
  - One shared Family legend (union of both genes' families) instead of two
    separate ones. Genus/Species legends stay separate per gene (their
    union is *larger* than either single-gene panel, so sharing wouldn't
    help there).
  - Sample x-tick labels have their redundant "-rbcl"/"-trnl" suffix
    stripped (each panel is already single-gene).
  - NA is always stacked last (topmost segment) in every bar, instead of
    wherever it happens to sort alphabetically among taxon names.
  - Legends are drawn in their own dedicated GridSpec rows rather than
    per-axes bbox_to_anchor legends.

Run from this directory (needs pandas + matplotlib, e.g. the
"metabarcoding" conda env on this machine):
    python make_combined_RRA_stacked_impala.py
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


def load_gene(csv_path):
    df = pd.read_csv(csv_path).fillna("NA")
    sample_columns = [c for c in df.columns if c not in ['ASV_ID', 'Family', 'Genus', 'Species', 'Life form']]
    df['Genus_Species'] = df.apply(
        lambda x: f"{x['Genus']} {x['Species']}" if x['Species'] != "NA" else f"{x['Genus']}", axis=1
    )
    return df, sample_columns


def strip_gene_suffix(sample_name):
    return sample_name.split('-')[0]


def per_sample_percent(df, sample_columns, level):
    grouped = df.groupby(level)[sample_columns].sum().T
    percent = grouped.div(grouped.sum(axis=1), axis=0) * 100
    percent.index = [strip_gene_suffix(s) for s in percent.index]
    if "NA" in percent.columns:
        ordered = [c for c in percent.columns if c != "NA"] + ["NA"]
        percent = percent[ordered]
    return percent


def is_gray(rgb, tol=0.1):
    return max(rgb) - min(rgb) < tol


def build_palette():
    all_colors = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
    return [c for c in all_colors if not is_gray(c)]


def build_shared_family_color_map(all_families, palette):
    families_sorted = sorted(f for f in all_families if f != "NA")
    color_map = {fam: palette[i % len(palette)] for i, fam in enumerate(families_sorted)}
    color_map["NA"] = "#B0B0B0"
    return color_map


def genus_of(genus_species_label):
    return genus_species_label.split()[0] if genus_species_label != "NA" else "NA"


def build_genus_color_map(genus_species_labels, palette):
    """Colour is assigned per GENUS (species under the same genus share a
    colour), so the legend only needs one entry per genus, not per
    genus+species combination.

    Genera are sorted alphabetically, which is also the order bars stack
    them in -- so sequential palette indices would land on alphabetically
    (and therefore stack-) adjacent genera. tab20/tab20b/tab20c each group
    similar shades at consecutive indices (tab20 pairs dark/light of the
    same hue back-to-back), so a plain `palette[i % n]` would frequently
    put two similar-looking shades right next to each other in a stacked
    bar. Using a stride coprime to the pair/quad grouping spreads
    sequential assignments across different hue families instead.
    """
    genera = sorted(set(genus_of(label) for label in genus_species_labels) - {"NA"})
    n = len(palette)
    stride = 7  # coprime to tab20's pairs (2) and tab20b/c's quads (4)
    color_map = {genus: palette[(i * stride) % n] for i, genus in enumerate(genera)}
    color_map["NA"] = "#B0B0B0"
    return color_map


def resolve_column_colors(columns, genus_color_map):
    """Map each Genus_Species column to its genus's colour, for passing to
    draw_stacked_panel (which expects one colour per column)."""
    return {c: genus_color_map[genus_of(c)] for c in columns}


def draw_stacked_panel(ax, percent_df, color_map, fontsize):
    cols = percent_df.columns
    percent_df.plot(kind='bar', stacked=True, ax=ax, color=[color_map[c] for c in cols],
                     legend=False, width=0.85)
    ax.set_ylabel('% Relative Read Abundance', fontsize=BODY_FONT_SIZE)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.tick_params(axis='y', labelsize=BODY_FONT_SIZE)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center', fontsize=fontsize)
    ax.set_xlabel('')


def draw_legend_row(ax, labels, color_map, ncol, fontsize, title=None):
    ax.axis('off')
    handles = [mpatches.Patch(color=color_map[label], label=str(label)) for label in labels]
    # mode='expand' forces the legend to fit exactly within this axes' width
    # instead of overflowing past the figure edge when there are many/long entries.
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0, 0, 1, 1), mode='expand',
              ncol=ncol, fontsize=fontsize, title=title, title_fontsize=fontsize,
              frameon=False, borderaxespad=0, labelspacing=0.3, columnspacing=0.8, handletextpad=0.4)


def legend_height(n_entries, ncol, row_pitch):
    nrows = -(-n_entries // ncol)  # ceil division
    return nrows * row_pitch + 0.15


def measure_text_width_in(text, fontsize):
    """Actual rendered width of `text` at `fontsize`, in inches -- measuring
    real glyph metrics instead of guessing avoids legend columns that are
    too narrow for the label and end up overlapping."""
    tmp_fig = plt.figure()
    renderer = tmp_fig.canvas.get_renderer()
    artist = tmp_fig.text(0, 0, text, fontsize=fontsize)
    width_in = artist.get_window_extent(renderer=renderer).width / tmp_fig.dpi
    plt.close(tmp_fig)
    return width_in


def compute_ncol(labels, fontsize, available_width_in, patch_width_in=0.22, gap_in=0.15):
    longest = max((str(l) for l in labels), key=len)
    entry_width = patch_width_in + measure_text_width_in(longest, fontsize) + gap_in
    return max(1, int(available_width_in // entry_width))


def main():
    rbcl_df, rbcl_cols = load_gene(HERE / "I_rbcl_revised" / "I_rbcl_diet_list_filtered_final.csv")
    trnl_df, trnl_cols = load_gene(HERE / "I_trnl_revised" / "I_trnl_diet_list_filtered_final.csv")

    rbcl_family_pct = per_sample_percent(rbcl_df, rbcl_cols, 'Family')
    trnl_family_pct = per_sample_percent(trnl_df, trnl_cols, 'Family')
    rbcl_genus_pct = per_sample_percent(rbcl_df, rbcl_cols, 'Genus_Species')
    trnl_genus_pct = per_sample_percent(trnl_df, trnl_cols, 'Genus_Species')

    palette = build_palette()
    shared_family_colors = build_shared_family_color_map(
        set(rbcl_family_pct.columns) | set(trnl_family_pct.columns), palette)
    shared_genus_labels = sorted(set(rbcl_genus_pct.columns) | set(trnl_genus_pct.columns))
    shared_genus_colors = build_genus_color_map(shared_genus_labels, palette)
    rbcl_genus_col_colors = resolve_column_colors(rbcl_genus_pct.columns, shared_genus_colors)
    trnl_genus_col_colors = resolve_column_colors(trnl_genus_pct.columns, shared_genus_colors)

    n_family = len(shared_family_colors)
    legend_genera = sorted(g for g in shared_genus_colors if g != "NA") + ["NA"]
    n_genus = len(legend_genera)

    # Column counts are computed from actually-measured text width, not guessed --
    # a guessed ncol that's too high makes long taxon names overlap each other.
    avail_w = FIG_WIDTH_IN * 0.95  # small safety margin below the nominal figure width
    # Both legends use the same, smallest allowed (Wiley floor) size for consistency.
    legend_fontsize = 6

    family_ncol = compute_ncol(shared_family_colors.keys(), legend_fontsize, avail_w)
    genus_ncol = compute_ncol(legend_genera, legend_fontsize, avail_w)

    row_pitch = 0.12
    family_legend_h = legend_height(n_family, ncol=family_ncol, row_pitch=row_pitch)
    genus_legend_h = legend_height(n_genus, ncol=genus_ncol, row_pitch=row_pitch)

    bar_row_h = 1.9

    # Each gap is set independently (rather than one shared GridSpec hspace)
    # since they need different amounts of room: the bar-row gap has to clear
    # the 90-degree-rotated sample labels hanging below each bar panel, while
    # the legend gaps can be much tighter.
    title_gap_in = 0.5
    bar_gap_in = 0.55       # between the two rows of bar panels (room for rotated x-tick labels + letter)
    section_gap_in = 0.35   # between the last bar row and the Family legend
    inter_legend_gap_in = 0.08
    bottom_gap_in = 0.05
    col_gap_in = 0.55
    left_margin_in = 0.4   # room for panel A/C's y-axis label, which has nothing to its left otherwise
    right_margin_in = 0.1  # small safety margin on the right edge

    fig_height = (title_gap_in + 2 * bar_row_h + bar_gap_in + section_gap_in
                  + family_legend_h + inter_legend_gap_in + genus_legend_h + bottom_gap_in)

    fig = plt.figure(figsize=(FIG_WIDTH_IN, fig_height))

    col_width_in = (FIG_WIDTH_IN - left_margin_in - right_margin_in - col_gap_in) / 2
    col0_x = left_margin_in
    col1_x = col0_x + col_width_in + col_gap_in

    def rect(x_in, y_top_in, w_in, h_in):
        return [x_in / FIG_WIDTH_IN, 1 - (y_top_in + h_in) / fig_height,
                w_in / FIG_WIDTH_IN, h_in / fig_height]

    y_row0 = title_gap_in
    y_row1 = y_row0 + bar_row_h + bar_gap_in
    y_family_legend = y_row1 + bar_row_h + section_gap_in
    y_genus_legend = y_family_legend + family_legend_h + inter_legend_gap_in

    ax_a = fig.add_axes(rect(col0_x, y_row0, col_width_in, bar_row_h))
    ax_b = fig.add_axes(rect(col1_x, y_row0, col_width_in, bar_row_h))
    ax_c = fig.add_axes(rect(col0_x, y_row1, col_width_in, bar_row_h))
    ax_d = fig.add_axes(rect(col1_x, y_row1, col_width_in, bar_row_h))

    draw_stacked_panel(ax_a, rbcl_family_pct, shared_family_colors, BODY_FONT_SIZE)
    draw_stacked_panel(ax_b, rbcl_genus_pct, rbcl_genus_col_colors, BODY_FONT_SIZE)
    draw_stacked_panel(ax_c, trnl_family_pct, shared_family_colors, BODY_FONT_SIZE)
    draw_stacked_panel(ax_d, trnl_genus_pct, trnl_genus_col_colors, BODY_FONT_SIZE)

    # Per-panel sub-titles removed -- the manuscript caption already states what
    # each panel shows (A/B = rbcL family/genus-species, C/D = trnL), so the
    # letters alone are enough and this saves the space/clashing the titles caused.
    for ax, letter in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        ax.text(-0.08, 1.05, letter, transform=ax.transAxes, fontsize=HEADING_FONT_SIZE, fontweight='bold')

    fig.suptitle("Figure A4", fontsize=HEADING_FONT_SIZE, y=0.985)

    # Just two legends now: one shared Family legend, one shared Genus legend
    # (species dropped from the legend -- colour is per-genus, and species
    # detail is already shown by the axis labels on the other RRA figures).
    ax_leg_family = fig.add_axes(rect(0.0, y_family_legend, FIG_WIDTH_IN, family_legend_h))
    ax_leg_genus = fig.add_axes(rect(0.0, y_genus_legend, FIG_WIDTH_IN, genus_legend_h))

    draw_legend_row(ax_leg_family, sorted(shared_family_colors), shared_family_colors,
                     ncol=family_ncol, fontsize=legend_fontsize, title="Family")
    draw_legend_row(ax_leg_genus, legend_genera, shared_genus_colors,
                     ncol=genus_ncol, fontsize=legend_fontsize, title="Genus")

    out_path = HERE / "Figure_A4_composite.png"
    fig.savefig(out_path, dpi=TARGET_DPI, bbox_inches='tight')

    with Image.open(out_path) as im:
        w_mm, h_mm = im.width / TARGET_DPI * 25.4, im.height / TARGET_DPI * 25.4
        size_kb = out_path.stat().st_size / 1024
        print(f"Saved {out_path}")
        print(f"  {im.width}x{im.height}px @ {TARGET_DPI} DPI -> {w_mm:.1f} x {h_mm:.1f} mm, {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
