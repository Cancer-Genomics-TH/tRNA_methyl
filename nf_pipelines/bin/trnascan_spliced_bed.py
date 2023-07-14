#!/usr/bin/env python3

"""
parse trnascan bed getting entries with 2 exons (column 10)
output: bed file for extracting spliced tRNAs

"""


import executor
import sys
import polars as pl

GENOME_FASTA = "chm13_2.0.fa"
OUTPUT_SPLD_BED = "chr_trnascan_spliced_trnas.bed"

def spliced_genes(input_bed):
    """
    takes the original tRNAScan-SE BED file
    outputs BED file with spliced tRNAs only
    """
    df = pl.read_csv(input_bed, separator="\t", has_header=False, infer_schema_length=10000).filter(pl.col("column_10") == 2)
    #debug
    #print(df)
    df.write_csv(OUTPUT_SPLD_BED, separator="\t", has_header=False)

def extract_spliced_trnas(OUTPUT_SPLD_BED):
    """
    takes the tRNAScan-SE BED file with spliced entries only
    outputs fasta 
    """
    
    command_1 = f"""bedtools getfasta -fi {GENOME_FASTA} -bed {OUTPUT_SPLD_BED} -fo {OUTPUT_SPLD_BED}.fa"""
    executor.execute(command_1, capture=True) 
    

if __name__ == "__main__":
    input_bed = sys.argv[1]
    spliced_genes(input_bed)
    extract_spliced_trnas(OUTPUT_SPLD_BED)
    

# Path: nf_pipelines/bin/trnascan_spliced_bed.py
