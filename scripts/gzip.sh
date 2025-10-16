#!/usr/bin/env bash
# Compress all .fastq files in the directory into .gz files : used on Illumina raw reads files to save space

INPUT_DIR="/data/projet3/Input/reads_paired_end/reads"    # dossier contenant les fichiers fastq
OUTPUT_DIR="/data/projet3/Input/reads_paired_end/compressed_reads.gz"  # dossier cible

for f in "$INPUT_DIR"/*.fastq; do   # loop through all .fastq files
    gzip -c "$f" > "$OUTPUT_DIR/$(basename "$f").gz"  # compress and save as .gz
    # -c : write to standard output and keep original file unchanged
    echo "Compressed: $f -> $OUTPUT_DIR/$(basename "$f").gz" 
done
