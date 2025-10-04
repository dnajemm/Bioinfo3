#Compress all .fastq files in the  directory into .zip files : used on Illumina raw reads
for f in *.fastq; do
    gzip -k "$f"
done
#-k option keeps the original .fastq files

