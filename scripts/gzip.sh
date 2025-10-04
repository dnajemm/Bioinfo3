# Compress all .fastq files in the directory into .gz files : used on Illumina raw reads
#!/usr/bin/env bash

INPUT_DIR="/data/projet3/reads_paired_end/reads"    # dossier contenant les fichiers fastq
OUTPUT_DIR="/data/projet3/reads_paired_end/compressed_reads.gz"  # dossier cible

for f in "$INPUT_DIR"/*.fastq; do   # loop through all .fastq files
    gzip -c "$f" > "$OUTPUT_DIR/$(basename "$f").gz"  # compress and save as .gz
    echo "Compressed: $f -> $OUTPUT_DIR/$(basename "$f").gz" 
done
