#!/bin/bash

OUT_DIR=/data/projet3/Input/reads_paired_end/reads  #output directory

FASTQ_DUMP=/data/projet3/sra/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump  #path to fastq-dump tool


while read -r SRA_ID; do
    $FASTQ_DUMP --split-files $SRA_ID -O $OUT_DIR  #--split-files to get paired-end reads
done < /data/projet3/Input/reads_paired_end/sra_list.txt   #input file with SRA IDs