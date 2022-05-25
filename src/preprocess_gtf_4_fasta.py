#!/usr/bin/env python3

"""
script for transforming a GTF with tRNA genes to a GTF with a gene_name in the second column
this fixes the problem:
when extracting fasta sequences using bedtools based on annotation

the name was not the gene name but value in the second column of the GTF

usage:

preprocess_gtf_4_fasta.py input.gtf > output.gtf

"""

import sys

input_gtf_fn = sys.argv[1]

with open(input_gtf_fn, encoding="utf-8") as input_fh:
    for line in input_fh:
        if line[0] == "#":
            pass
        else:
            line = line.strip()
            sl = line.split("\t")
            # print(sl)
            col_nine_split = sl[8].split(" ")
            gene_name = col_nine_split[5]
            gene_name = gene_name[1:-2]
            ens_id = col_nine_split[1]
            out_str_1 = "\t".join(sl[:2])
            out_str_2 = "\t".join(sl[3:])
            out_str_1 = f"{out_str_1}\t{gene_name}"
            out_str = f"{out_str_1}\t{out_str_2}"
            print(out_str)
