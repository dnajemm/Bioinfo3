#!/bin/bash

INPUT_DIR="/data/projet3/Input/reads_paired_end/compressed_reads.gz" #input directory
OUTPUT_DIR="/data/projet3/Output/fastqc_results" #output directory

for file in "$INPUT_DIR"/*.fastq.gz ; do  #for each file in the input directory
        fastqc -o "$OUTPUT_DIR" "$file" #run fastqc and output to the output directory
done
