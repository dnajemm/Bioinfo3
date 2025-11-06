#!/usr/bin/env python3
"""The script counts the total number of SNPs across all VCF files in a specified folder and provides a breakdown of SNP counts per strain.
It outputs the total SNP count, a per-strain SNP count, and saves the results as a CSV file and a PNG image of a table.
Usage: python count_number_SNPs.py <folder_with_vcf_files>
"""
import sys, os, gzip

folder = sys.argv[1]
total_variants = 0
strain_counts = {} 
# For each vcf file in the folder given as argument
for f in os.listdir(folder): 
    if f.endswith('.vcf') or f.endswith('.vcf.gz'):
        path = os.path.join(folder, f)
        op = gzip.open if f.endswith('.gz') else open
        with op(path, 'rt') as fh:
            # SUM the number of variants of this file
            n = sum(1 for l in fh if not l.startswith('#'))
            total_variants += n
            strain_name = f.replace('.vcf.gz', '').replace('.vcf', '').replace('_filtered', '')
            strain_counts[strain_name] = n

# Return the total number of all the variants found for all strains present : 
print(total_variants)
# Return the number of variants per strain :
print("\nSNP count per strain:")
for strain, count in sorted(strain_counts.items()):
    print(f"{strain}: {count}")

# Export results as CSV + PNG :
import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame(list(strain_counts.items()), columns=['Strain', 'SNP_Count'])

# save as PNG table
fig, ax = plt.subplots(figsize=(8, len(strain_counts)*0.3 + 1))
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5) 
# Style the header row to add a colored background and bold text
for (row, col), cell in table.get_celld().items():
    if row == 0:  # header row
        cell.set_facecolor("#bcd4e6")  # light blue background
        cell.set_text_props(weight='bold', color='black')
        cell.set_fontsize(11)
        cell.set_height(cell.get_height() * 1.1)

# Save both files in the output directory
out_dir = '/data/projet3/Output/results_numbers_SNP_het'
os.makedirs(out_dir, exist_ok=True)

csv_path = os.path.join(out_dir,'SNP_counts_per_strain.csv')
df.to_csv(csv_path, index=False)

png_path = os.path.join(out_dir,'SNP_counts_per_strain.png')
plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.close()
