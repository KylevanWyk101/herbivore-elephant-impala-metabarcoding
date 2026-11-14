"""
Shannon diversity boxplot+stripplot (Figure 1), resized to fit the Wiley
Ecology and Evolution appendix/figure spec (<=18x22cm, 300-600 DPI).

The original figsize=(12,6) in significance_test_mann_whitney_only.py
produces a 300x147mm page -- fine resolution-wise (this is a native
matplotlib vector PDF, no PowerPoint step, no rasterization loss involved),
just physically too wide. The plotting logic below (boxplot, stripplot,
legend, p-value annotations) is copied verbatim from that script; only the
figure size and the added "Figure 1" title are new. That script, its CSV,
and its own outputs are untouched by this one.

Run from this directory (needs pandas, seaborn, scipy, matplotlib):
    python make_figure1_composite.py
"""
from pathlib import Path

from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent

TARGET_DPI = 1000  # top of Wiley's 600-1000 DPI "line art" preferred range (electronic_artwork_guidelines.pdf)
FIG_WIDTH_IN = 170 / 25.4  # under Wiley's 180mm cap
FIG_HEIGHT_IN = FIG_WIDTH_IN * (6 / 12)  # same aspect ratio as the original figsize=(12,6)

df = pd.read_csv(HERE / "combined_shannon_index_genus.csv")
df['Group'] = df['Group'].str.strip().str.capitalize()
df['Locus'] = df['Locus'].str.strip().str.lower()

locus_labels = {
    'rbcl': r'$\it{rbc}$L',
    'trnl': r'$\it{trn}$L'
}
df['Locus_display'] = df['Locus'].map(locus_labels).fillna(df['Locus'])

unique_loci = df['Locus'].unique()
results = []
for locus in unique_loci:
    sub_df = df[df['Locus'] == locus]
    elephant = sub_df[sub_df['Group'] == 'Elephant']['ShannonIndex']
    impala = sub_df[sub_df['Group'] == 'Impala']['ShannonIndex']
    stat, p_val = mannwhitneyu(elephant, impala, alternative='two-sided')
    results.append({'Locus': locus, 'Test': "Mann-Whitney U", 'Statistic': stat, 'p-value': p_val})

BODY_FONT_SIZE = 7   # Wiley requires figure text sized 6-10pt at final print size
HEADING_FONT_SIZE = 9

# font_scale is relative to seaborn's ~10pt base, so 0.7 lands text close to
# BODY_FONT_SIZE -- the original font_scale=1.1 was sized for the original
# 12x6in canvas and became way oversized once the canvas shrunk to fit Wiley's
# width cap (absolute point sizes don't shrink automatically with figsize).
sns.set(style="whitegrid", font_scale=0.7)
# sns.set() bundles its own grid/spine widths into the "whitegrid" style, so
# these have to be set after it (an earlier assignment would just get
# overwritten) to actually take effect and match the thinner lines used in
# the other figures in this project.
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['grid.linewidth'] = 0.4
fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

ax = sns.boxplot(data=df, x='Locus_display', y='ShannonIndex', hue='Group',
                  palette='pastel', showfliers=False, width=0.5, linewidth=0.8)

sns.stripplot(data=df, x='Locus_display', y='ShannonIndex', hue='Group',
              dodge=True, marker='o', alpha=0.7, size=4,
              linewidth=0.4, palette='dark', edgecolor='gray')

handles, labels = ax.get_legend_handles_labels()
# Pinned to upper-left (matching the original figure's layout) rather than
# left to 'best' auto-placement, which moved it around and needed the axis
# stretched off-centre to make room -- that broke the plot's symmetry.
legend = ax.legend(handles[:2], labels[:2], title="Group", loc='upper left',
                    fontsize=BODY_FONT_SIZE, title_fontsize=BODY_FONT_SIZE,
                    markerscale=0.7, handletextpad=0.4, borderpad=0.4, labelspacing=0.3)
legend.get_frame().set_linewidth(0.5)
for patch in legend.legend_handles:
    patch.set_linewidth(0.4)

for i, locus in enumerate(unique_loci):
    p_val = results[i]['p-value']
    label = f"$\\mathit{{p}}$={p_val:.3f}" if p_val >= 0.05 else "$\\mathit{{p}}$<0.05"

    group_max = df[df['Locus'] == locus]['ShannonIndex'].max()
    y_offset = 0.05 * (df['ShannonIndex'].max() - df['ShannonIndex'].min())
    y_pos = group_max + y_offset

    ax.text(i, y_pos, label, ha='center', va='bottom', fontsize=BODY_FONT_SIZE, color='black')

plt.ylabel("Shannon Diversity Index")
plt.xlabel("Locus")
fig.suptitle("Figure 1", fontsize=HEADING_FONT_SIZE, y=1.02)
plt.tight_layout()

# Render at 2x TARGET_DPI internally, then downsample with a high-quality
# filter to the pixel size TARGET_DPI actually implies -- on top of already
# raising TARGET_DPI itself to 1000 (the top of the line-art range), this
# gets the sharpest result the guidelines allow. Small circular markers
# (mostly antialiased edge, little fill) benefit from this the most.
SUPERSAMPLE = 2
out_path = HERE / "Figure_1_composite.png"
plt.savefig(out_path, format="png", dpi=TARGET_DPI * SUPERSAMPLE, bbox_inches="tight")

with Image.open(out_path) as im:
    target_size = (round(im.width / SUPERSAMPLE), round(im.height / SUPERSAMPLE))
    im = im.resize(target_size, Image.LANCZOS)
    im.save(out_path, dpi=(TARGET_DPI, TARGET_DPI), optimize=True)
    w_mm, h_mm = im.width / TARGET_DPI * 25.4, im.height / TARGET_DPI * 25.4
    size_kb = out_path.stat().st_size / 1024
    print(f"Saved {out_path}")
    print(f"  {im.width}x{im.height}px @ {TARGET_DPI} DPI -> {w_mm:.1f} x {h_mm:.1f} mm, {size_kb:.0f} KB")

# Also save as a vector PDF -- the guidelines list EPS/PDF as the *preferred*
# format for line art (PNG is only "acceptable"). Unlike the PNG, this has no
# pixel grid for the actual lines/text/markers, so it can't look fuzzy at any
# zoom level -- a more fundamental fix than raising raster DPI further.
pdf_path = HERE / "Figure_1_composite.pdf"
plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
print(f"Saved {pdf_path} (vector, {pdf_path.stat().st_size / 1024:.0f} KB)")
