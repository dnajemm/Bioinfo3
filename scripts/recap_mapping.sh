#!/bin/bash
# Script to generate a summary CSV file with mapping statistics for BAM files.
# Takes as input a folder containing BAM files.
# Outputs a CSV file named recap_mapping.csv with columns:
# Fichier,Taux_mapping(%),Profondeur_moyenne

# Define input and output directories
BAM_DIR="/data/projet3/Output/results_mapping"
OUT_CSV="$BAM_DIR/recap_mapping.csv"

echo "Fichier,Taux_mapping(%),Profondeur_moyenne" > "$OUT_CSV"

# Loop through all BAM files in the specified folder 
for bam in "$BAM_DIR"/*.bam
do
    echo "Traitement de $bam ..."
    
    #  total number of reads
    total=$(samtools flagstat "$bam" | head -n 1 | awk '{print $1}')
    
    # number of mapped reads (line "mapped (")) 
    mapped=$(samtools flagstat "$bam" | grep " mapped (" | head -n 1 | awk '{print $1}')
    
    # % of mapped reads
    rate=$(echo "scale=2; 100 * $mapped / $total" | bc)
    
    # Average depth calculation
    mean_depth=$(samtools depth "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    
    # Write to CSV file
    echo "$(basename "$bam"),$rate,$mean_depth" >> "$OUT_CSV"
done


