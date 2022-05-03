#!/usr/bin/env pypy3

"""
remove "chr" prefix from chromosome names
rename "chrM" to ENSEMBL compatible "MT"

using pypy3 for speed

usage:

fix_chr_names.py > GRCh38.primary_assembly.genome.names_fix.fa

"""

from pyfaidx import Fasta

genome_fn = "GRCh38.primary_assembly.genome.fa"

gencode_genome = Fasta(genome_fn)
for seq_id in gencode_genome.keys():
    if seq_id[:3] == "chr":
        if seq_id == "chrM":
            new_seq_id = "MT"
        else:
            new_seq_id = seq_id[3:]
    else:
        new_seq_id = seq_id
    fasta_id = f">{new_seq_id}"
    print(fasta_id)
    for line in gencode_genome[seq_id]:
        print(line)
