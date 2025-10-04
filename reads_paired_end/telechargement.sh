#!/bin/bash

OUT_DIR=/data/projet3/reads_paired_end/

FASTQ_DUMP=/data/projet3/sra/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump


while read -r SRA_ID; do
    echo "Téléchargement de $SRA_ID..."
    $FASTQ_DUMP --split-files $SRA_ID -O $OUT_DIR
done < /data/projet3/reads_paired_end/sra_list.txt