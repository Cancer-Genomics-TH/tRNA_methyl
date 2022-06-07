#!/usr/bin/env python3

"""
convert tab delim file from gtRNAdb to a stockholm format aligment

GGGGGTATAGCTCAGT.GGTA.GAGCGCGTGCTTAGCATGCACGAG...........GtCCTGGGTTCGATCCCCAGTACCTCCA... tRNA-Ala-AGC-1-1

input file: 
    meta_data/tRNA_numbering/gtrna_align.isotypes.all.2cols

usage:
    ./tab_gtrna_to_stockholm.py gtrna_align.isotypes.all.2cols > gtrna_align.isotypes.all.sto

"""

import sys

align_header = "# STOCKHOLM 1.0"


input_fn = sys.argv[1]

with open(input_fn) as f:
    old_id_prefix = ""
    counter = 0
    for line in f:
        line = line.strip()
        seq, seq_id = line.split()
        split_id = seq_id.split("-")
        id_prefix = "-".join(split_id[:3])

        if id_prefix != old_id_prefix:
            if counter != 0:
                print("//")
            print(align_header)
            align_group_id = f"#=GF ID {id_prefix}"
            print(align_group_id)
        else:
            pass
        print(f"{seq_id}\t{seq}")
        old_id_prefix = id_prefix
        counter += 1

print("//")
