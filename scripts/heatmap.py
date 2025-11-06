"""
The script cleans genomic data, converts LOH/HET to numeric values,
divides each chromosome into 1 kb bins to compare all strains on the same grid.
It plots one heatmap per chromosome, and a combined genome-wide heatmap
with colored bars for clusters (Y-axis) and chromosomes (top).
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# Paths

INPUT = "/data/projet3/Output/results_windows/all_samples.het.5kb_LOH_recalc.tsv"
OUT = Path("/data/projet3/Output/results_heatmap")
CLUSTERS_XLSX = "/data/projet3/Input/tableauSouches.xlsx"
OUT.mkdir(parents=True, exist_ok=True)

# Load + clean data

df = pd.read_csv(INPUT, sep="\t", usecols=["SAMPLE", "CHROM", "START", "END", "STATUS"])
df["START"] = pd.to_numeric(df["START"], errors="coerce")
df["END"] = pd.to_numeric(df["END"], errors="coerce")
df = df.dropna(subset=["START", "END"])
df = df[df["END"] > df["START"]]

# Load clusters
clusters = pd.read_excel(CLUSTERS_XLSX).rename(columns={"Isolate": "SAMPLE", "Cluster": "CLUSTER"})
clusters["SAMPLE"] = clusters["SAMPLE"].astype(str).str.strip()
clusters["CLUSTER"] = clusters["CLUSTER"].astype(str).str.strip()
df["SAMPLE"] = df["SAMPLE"].astype(str).str.strip()

# Merge cluster info
df = df.merge(clusters[["SAMPLE", "CLUSTER"]], on="SAMPLE", how="left")

# Map LOH/HET to numeric
df["LOH"] = df["STATUS"].astype(str).str.strip().str.upper().map({"HET": 0.0, "LOH": 1.0})

# Create 1 kb bins across chromosomes
STEP = 1000
rows = []
for _, row in df.iterrows():
    start, end = int(row["START"]), int(row["END"])
    if end <= start:
        continue
    for p in range((start // STEP) * STEP, ((end - 1) // STEP) * STEP + STEP, STEP):
        rows.append([row["SAMPLE"], row["CHROM"], p, row["LOH"], row["CLUSTER"]])

bins = pd.DataFrame(rows, columns=["SAMPLE", "CHROM", "POS", "LOH", "CLUSTER"])
agg = bins.groupby(["SAMPLE", "CHROM", "POS"], as_index=False)["LOH"].max()
agg = agg.merge(bins[["SAMPLE", "CLUSTER"]].drop_duplicates(), on="SAMPLE", how="left")

# Plot per-chromosome heatmaps

for chrom, sub in agg.groupby("CHROM"):
    mat = sub.pivot(index="SAMPLE", columns="POS", values="LOH")
    if mat.empty:
        continue
    mat = mat.reindex(sorted(mat.columns), axis=1)

    order = (
        sub[["SAMPLE", "CLUSTER"]]
        .drop_duplicates()
        .sort_values(["CLUSTER", "SAMPLE"], na_position="last"))
    mat = mat.loc[order["SAMPLE"]]
    # Prepare data for plotting
    A = mat.to_numpy(float)
    pos = mat.columns.to_numpy(int)
    x_edges = np.r_[pos, pos[-1] + STEP]
    y_edges = np.arange(mat.shape[0] + 1)
    # Plot heatmap
    plt.figure(figsize=(30, 10))
    plt.pcolormesh(x_edges, y_edges, A, shading="flat", cmap="bwr", vmin=0, vmax=1)
    plt.title(f"LOH heatmap — {chrom}")
    plt.xlabel("Genomic position (bp)")
    plt.ylabel("") 
    plt.yticks(np.arange(mat.shape[0]) + 0.5, mat.index, fontsize=12)  # sample names
    plt.tight_layout()
    plt.savefig(OUT / f"LOH_heatmap_{chrom}.png", dpi=200)
    plt.close()

# Combine all chromosomes side-by-side

chrom_order = sorted(agg["CHROM"].unique(), key=str)
combined, bounds, offset = [], [], 0

for chrom in chrom_order:
    sub = agg[agg["CHROM"] == chrom]
    mat_c = sub.pivot(index="SAMPLE", columns="POS", values="LOH")
    if mat_c.empty:
        continue
    mat_c = mat_c.reindex(sorted(mat_c.columns), axis=1)
    ncols = mat_c.shape[1]
    new_cols = np.arange(offset, offset + ncols * STEP, STEP)
    mat_c.columns = new_cols
    bounds.append((offset, offset + ncols * STEP, chrom))
    offset += ncols * STEP
    combined.append(mat_c)

mat_all = pd.concat(combined, axis=1)

# Order by cluster
order = (
    agg[["SAMPLE", "CLUSTER"]]
    .drop_duplicates()
    .sort_values(["CLUSTER", "SAMPLE"], na_position="last"))
mat_all = mat_all.loc[order["SAMPLE"]]
cluster_map = dict(zip(order["SAMPLE"], order["CLUSTER"]))

# Prepare cluster color mapping
cluster_levels = [c for c in order["CLUSTER"].unique() if pd.notna(c)]
cluster_palette = plt.cm.Set2.colors
cluster_colors = {c: cluster_palette[i % len(cluster_palette)] for i, c in enumerate(cluster_levels)}
cluster_colors[np.nan] = (0.85, 0.85, 0.85)

# Plot combined genome-wide heatmap
A = mat_all.to_numpy(float)
pos = np.array(sorted(mat_all.columns))
x_edges = np.r_[pos, pos[-1] + STEP]
y_edges = np.arange(mat_all.shape[0] + 1)
ylabels = list(mat_all.index)   # only strain names
fig, ax = plt.subplots(figsize=(40, 8))
pcm = ax.pcolormesh(x_edges, y_edges, A, shading="flat", cmap="bwr", vmin=0, vmax=1)
ax.set_facecolor("white")
ax.set_title("Genome-wide LOH heatmap (all chromosomes combined)")
ax.set_xlabel("Genomic position (concatenated chromosomes)")
ax.set_ylabel("")
ax.set_yticks(np.arange(mat_all.shape[0]) + 0.5)
ax.set_yticklabels(ylabels, fontsize=8)

# Colorbar for LOH/HET
cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label("LOH (1) / HET (0)", fontsize=16, fontweight="bold",labelpad=25)  
cbar.ax.tick_params(labelsize=16)  

# Add cluster color bar on the left (rectangles)
nrows = mat_all.shape[0]
row_cluster = [cluster_map.get(s, np.nan) for s in mat_all.index]
row_colors = [cluster_colors.get(c, (0.85, 0.85, 0.85)) for c in row_cluster]
# Add inset axes for cluster colors
ax_left = ax.inset_axes([-0.025, 0, 0.02, 1])
ax_left.set_ylim(0, nrows)
ax_left.set_xlim(0, 1)
ax_left.axis("off")
# Draw rectangles for each row
for i, col in enumerate(row_colors):
    ax_left.add_patch(Rectangle((0, i), 1, 1, facecolor=col, edgecolor="none"))
# Label for clusters
ax_left.text(-0.2, nrows / 2, "Clusters", rotation=90, va="center", ha="right",
             fontsize=20, transform=ax_left.transData)

# Add chromosome color strip at the top
ax_top = ax.inset_axes([0, 1.02, 1, 0.08])
ax_top.set_xlim(ax.get_xlim())
ax_top.set_ylim(0, 1)
ax_top.axis("off")
# Draw rectangles for each chromosome
palette = plt.cm.tab20.colors
chr_colors = []
for i, (start, end, chrom) in enumerate(bounds):
    color = palette[i % 20]
    chr_colors.append((chrom, color))
    ax_top.add_patch(Rectangle((start, 0), end - start, 1, facecolor=color, edgecolor="none"))
    if i > 0:
        ax.axvline(start, color="white", linewidth=1)
# Label for chromosomes
ax_top.text(
    0.5 * sum(ax_top.get_xlim()),
    0.5,
    "Chromosomes",
    ha="center", va="center",
    fontsize=20)

# Chromosome legend handles
chr_handles = [
    Line2D([0], [0],
           marker="s", linestyle="",
           markerfacecolor=color, markeredgecolor="k",
           markersize=12, label=str(chrom))
    for chrom, color in chr_colors]
# Chromosome legend 
ax_top.legend(
    handles=chr_handles,
    title="Chromosomes",
    loc="upper left",
    bbox_to_anchor=(1.15, 0.30),   
    frameon=False,
    fontsize=16,
    title_fontsize=16,
    handlelength=1.5,
    handleheight=1.5,
    labelspacing=1.0)

# Clusters :
# sort clusters naturally
_raw = [c for c in order["CLUSTER"].unique() if pd.notna(c)]
def _natkey(x):
    s = str(x).strip()
    try:
        return (0, int(s))   # number first
    except ValueError:
        return (1, s.lower())  # then alphanum
cluster_levels_sorted = sorted(_raw, key=_natkey)

# handles for the legend
cluster_handles = [
    Line2D([0], [0],
           marker="s", linestyle="",
           markerfacecolor=cluster_colors[c], markeredgecolor="k",
           markersize=12, label=str(c))
    for c in cluster_levels_sorted]

# reserve space on the right for the legend
fig.subplots_adjust(right=0.86)

# legend on the right, properly aligned
leg = fig.legend(
    handles=cluster_handles,
    title="Clusters",
    loc="center left",
    bbox_to_anchor=(0.85, 0.50),
    frameon=False,
    borderaxespad=0.0,
    fontsize=14,           
    title_fontsize=16,    
    handlelength=1.5,     
    handleheight=1.5,     
    labelspacing=1.2)

# Adjust layout
fig.tight_layout()
fig.subplots_adjust(left=0.08, right=0.83)  
ax.set_title("Genome-wide LOH heatmap (all chromosomes combined)", fontsize=20, fontweight="bold")
ax.set_xlabel("Genomic position (concatenated chromosomes)", fontsize=20)
# Save figure
fig.savefig(OUT / "LOH_heatmap_combined_with_chr_cluster_legends.png")
plt.close(fig)
