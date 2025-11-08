#!/usr/bin/env python3
import glob

# Tous les fichiers 5 kb à traiter
files = glob.glob("/data/projet3/Output/results_windows/*.het.5kb.tsv")

# Seuil : LOH si N_HET < 3
threshold = 3

for input_file in files:

    output_file = input_file.replace(".het.5kb.tsv", ".het.5kb_LOH.tsv")

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for i, line in enumerate(fin):
            line = line.strip()

            # Entête -> on ajoute STATUS
            if i == 0:
                fout.write(line + "\tSTATUS\n")
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                # ligne invalide, on saute
                continue

            # 4e colonne = N_HET (entier)
            try:
                n_het = int(parts[3])
            except ValueError:
                # si valeur non numérique, par sécurité on met HET
                n_het = 999999

            status = "LOH" if n_het <= threshold else "HET"
            fout.write("\t".join(parts) + f"\t{status}\n")




#now for all 
awk 'BEGIN{OFS="\t"} NR==1{print $1,$2,$3,$4,$5,$6,"STATUS"; next} {print $1,$2,$3,$4,$5,$6,($5<=3?"LOH":"HET")}' \
/data/projet3/Output/results_windows/all_samples.het.5kb_LOH.tsv \
> /data/projet3/Output/results_windows/all_samples.het.5kb_LOH_recalc.tsv
# i just corrected with awk the last colonne 