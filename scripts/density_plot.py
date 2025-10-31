import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load all 5kb TSV files and concatenate them into a single DataFrame
files = sorted(glob.glob("/data/projet3/Output/results_windows/*.het.5kb.tsv"))
dfs = [pd.read_csv(f, sep=r"\s+|\t", engine="python", header=None,
                   names=["Chromosome","Start","End","HetSNPs"]) for f in files]
df = pd.concat(dfs, ignore_index=True)

# Clean and preprocess data trnansforming columns to numeric and dropping invalid rows
for col in ["Start","End","HetSNPs"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna().reset_index(drop=True)

# Calculate SNP density per kb
# Filter out invalid coordinates first
df = df[df["End"] > df["Start"]].copy()
win_lengths = df["End"] - df["Start"] + 1
df["SNP_density_per_kb"] = df["HetSNPs"] / win_lengths * 1000.0

# Remove negative densities
df = df[df["SNP_density_per_kb"] >= 0].copy()

# Create figure  
plt.figure(figsize=(10, 6))

# Create KDE plot 
sns.kdeplot(data=df, x="SNP_density_per_kb", linewidth=2.5)

# Add titles and labels
plt.title(f"Distribution of SNP density considering all {len(files)} strains", 
          fontsize=14)
plt.xlabel("SNPs per kb", fontsize=12)
plt.ylabel("Density", fontsize=12)

plt.tight_layout()
plt.savefig("/data/projet3/Output/result_density_plot/snp_density_distribution_all_strains.png")
plt.show()
