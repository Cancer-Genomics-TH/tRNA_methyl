#!/usr/bin/env pypy3

"""
remove "chr" prefix from chromosome names
rename "chrM" to ENSEMBL compatible "MT"

using pypy3 for speed

usage:

fix_chr_names.py > GRCh38.primary_assembly.genome.names_fix.fa

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
