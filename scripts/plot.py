import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Get all 5kb.tsv files
files = sorted(glob.glob("/data/projet3/Output/results_windows/*5kb.tsv"))

# Iterate over each file
for file in files:
    sample_name = os.path.basename(file).replace(".het.5kb.tsv", "")
    print(f"Processing {sample_name}...")

    df = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "HetSNPs"])
    df["Mid"] = (df["Start"] + df["End"]) / 2

    plt.figure(figsize=(30, 6), dpi=150)
    chroms = df["Chromosome"].unique()
    colors = plt.cm.tab20.colors

    # Trace every chromosome in color
    for i, chrom in enumerate(chroms):
        sub = df[df["Chromosome"] == chrom]     # filtering for the chromosome
        plt.plot(sub["Mid"], sub["HetSNPs"],    # plotting the data
                 color=colors[i % len(colors)], # cycle through colors if more than available
                 label=chrom, linewidth=1.2)    

    plt.xlabel("Position on chromosome (bp)")
    plt.ylabel("Number of heterozygous SNPs")
    plt.title(f"Heterozygosity profile — {sample_name}")
    plt.legend(title="Chromosome", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False) # legend outside plot ,bbox_to_anchor used to position it and loc to set location
    plt.grid(True, linestyle="--", alpha=0.5) #alternative grid style , alpha for transparency
    plt.tight_layout() # tight layout

    # saving the plot
    out_path = f"/data/projet3/Output/results_plot/plot_5Kb/{sample_name}_hetSNPs_all_chromosomes.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

print("Tous les graphiques 5kb ont été générés")

#For 10kb and 20kb, similar code blocks would follow with appropriate changes to file paths and names.