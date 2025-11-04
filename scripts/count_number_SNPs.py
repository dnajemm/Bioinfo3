#!/usr/bin/env python3
"""The script counts the total number of SNPs across all VCF files in a specified folder."""
import sys, os, gzip

folder = sys.argv[1]
total_variants = 0

# For each vcf file in the folder given as argument
for f in os.listdir(folder): 
    if f.endswith('.vcf') or f.endswith('.vcf.gz'):
        path = os.path.join(folder, f)
        op = gzip.open if f.endswith('.gz') else open
        with op(path, 'rt') as fh:
            # SUM the number of variants of this file
            total_variants += sum(1 for l in fh if not l.startswith('#')) 
# Return the total number of all the variants found for all strains present : 
print(total_variants)

