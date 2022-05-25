#!/usr/bin/env pypy3

"""
extracts Pseudo_tRNA sequences from a fasta

purpose: 
input for a check if sequence clustering of all tRNAs puts some of real tRNA sequences into contigs labeled as Pseudo_tRNA 

"""
import sys

from pyfaidx import Fasta

input_fasta_fn = sys.argv[1]

input_fasta = Fasta(input_fasta_fn)
for seq_id in input_fasta.keys():
    if seq_id[:11] == "Pseudo_tRNA":
        pass
    else:
        new_seq_id = seq_id
        fasta_id = f">{new_seq_id}"
        print(fasta_id)
        for line in input_fasta[seq_id]:
            print(line)
