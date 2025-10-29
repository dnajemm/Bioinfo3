# Script to plot heterozygosity profiles from TSV files , for 5kb, 10kb, and 20kb window sizes in a grid of subplots for each chromosome.
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import math
import numpy as np

# Get all 5kb.tsv files
# glob pattern to match files in the specified directory 
# sorted to ensure consistent order
files = sorted(glob.glob("/data/projet3/Output/results_windows/*5kb.tsv"))

# Iterate over each file
for file in files:
    sample_name = os.path.basename(file).replace(".het.5kb.tsv", "")
    df = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "HetSNPs"])
    df["Mid"] = ((df["Start"] + df["End"]) / 2)/1000 # Convert to kb for better readability

    # Get unique chromosomes
    chroms = df["Chromosome"].unique()
    n = len(chroms)
    ncols = 4
    nrows = math.ceil(n / ncols) # nbr of rows needed

    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(50, 10*nrows/ncols), dpi=150, sharey=True)
    axes = np.array(axes).reshape(-1) # flatten in case of single row
    colors = plt.cm.tab20.colors       # use a colormap with enough distinct colors

    # Plot each chromosome
    for i, chrom in enumerate(chroms):
        ax = axes[i]
        sub = df[df["Chromosome"] == chrom].sort_values("Mid")
        # Plot heterozygous SNPs vs position on chromosome 
        ax.plot(sub["Mid"], sub["HetSNPs"], 
                color=colors[i % len(colors)], linewidth=1.2) # cycle through colors if more chromosomes than colors
        ax.set_title(chrom) # set chromosome as title
        ax.grid(True, linestyle="--", alpha=0.5)
        if i // ncols == nrows - 1: # only bottom row
            ax.set_xlabel("Position on chromosome (kb)")
        if i % ncols == 0:
            ax.set_ylabel("Number of heterozygous SNPs")

    # Hide unused axes if the number of chromosomes is not a multiple of ncols
    for j in range(i+1, len(axes)):
        axes[j].axis("off")

    # Add overall title and adjust layout
    fig.suptitle(f"Heterozygosity profile — {sample_name}", y=1.02, fontsize=7)
    plt.tight_layout()

    # Save the figure
    out_path = f"/data/projet3/Output/results_plot/plot_5kb/{sample_name}_hetSNPs.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

#For 10kb and 20kb, similar code blocks would follow with appropriate changes to file paths and names.