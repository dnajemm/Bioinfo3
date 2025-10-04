#!/usr/bin/env bash

REF="/data/projet3/reference/brbr.fasta"           # reference genome
READS_DIR="/data/projet3/reads_paired_end/compressed_reads.gz"   # folder with *_1.fastq.gz and *_2.fastq.gz
OUT_DIR="/data/projet3/results_mapping"            # NEW output folder
THREADS=8

# Create the output folder (if it doesn’t exist)
mkdir -p "$OUT_DIR"
shopt -s nullglob   # ignore if no files match

for R1 in "$READS_DIR"/*_1.fastq.gz; do
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"               # pair file
  PREFIX="$(basename "$R1" _1.fastq.gz)"           # ex: YJS8112
  BAM="${OUT_DIR}/${PREFIX}.bam"

  echo ">> Mapping $PREFIX"
  bwa mem -t "$THREADS" -M \
    -R "@RG\tID:${PREFIX}\tSM:${PREFIX}" \
    "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "$BAM" -

  samtools index "$BAM"
done

echo "✅ Done. BAMs are in $OUT_DIR"
