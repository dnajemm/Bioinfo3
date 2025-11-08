#!/usr/bin/env python3
import glob

# Path to your results folder
files = glob.glob("/data/projet3/Output/results_windows/*.tsv")

for f in files:
    # Skip any merged ("all") files
    if "all" in f:
        continue

    # Choose the header based on the filename
    if "5kb_LOH" in f:
        header = "CHROM\tSTART\tEND\tN_HET\tSTATUS\n"
    elif "10kb" in f or "20kb" in f:
        header = "CHROM\tSTART\tEND\tN_HET\n"
    else:
        continue  # skip other file types

    # Read current content
    with open(f, "r") as fin:
        lines = fin.readlines()

    # Skip if header already exists
    if lines and lines[0].startswith("CHROM"):
        print(f"✔ Header already present: {f}")
        continue

    # Rewrite file with the header on top
    with open(f, "w") as fout:
        fout.write(header)
        fout.writelines(lines)
#i did a mistake on the second line so i corrected later on 
for f in /data/projet3/Output/results_windows/*.het.5kb_LOH.tsv; do
    awk 'BEGIN{OFS="\t"} 
         NR==1 {
             print "CHROM","START","END","N_HET","STATUS";  # clean header
             next
         } 
         {
             n=$4;
             status = (n <= 3 ? "LOH" : "HET");
             print $1,$2,$3,$4,status
         }' "$f" > "${f%.tsv}_fixed.tsv"
done
  
