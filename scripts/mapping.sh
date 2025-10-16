#!/usr/bin/env bash
# Script de mapping des lectures Illumina sur le génome de référence avec BWA et Samtools

REF="/data/projet3/Input/reference/brbr.fasta"            #reference genome
READS_DIR="/data/projet3/Input/reads_paired_end/compressed_reads.gz"   #folder with input reads
OUT_DIR="/data/projet3/Output/results_mapping"             # output folder
THREADS=8                                           # number of CPU cores used

mkdir -p "$OUT_DIR"                                # create output directory if it doesn't exist

for R1 in "$READS_DIR"/*_1.fastq.gz; do            # loop through all R1 files
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"               # corresponding R2 file
  PREFIX="$(basename "$R1" _1.fastq.gz)"           # sample name prefix
  BAM="${OUT_DIR}/${PREFIX}.bam"                   #  output BAM file

  bwa mem -t "$THREADS" -M \                       # align reads with BWA-MEM
    -R "@RG\tID:${PREFIX}\tSM:${PREFIX}" \         # read group information
    "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "$BAM" -        # sort and output BAM file

  samtools index "$BAM"                            # index the BAM file
done

