#!/usr/bin/env pypy3

"""
script to transform 2 lines of pre-parsed maf aligment to a fasta file
with a fasta names containing:

tRNA contig name and the number of counts
the sequence is from the NGS read match (gaps introduced during the mapping are removed )

uses pypy3 for speed

usage:

maf_2_fa_with_counts.py input.parsed.maf > output.fasta

"""


import string
import sys

FLOWCELL_PREFIX = "SND"


input_maf_fn = sys.argv[1]

sequence_frequency_dict = {}

with open(input_maf_fn, encoding="utf-8") as maf_fh:
    # old_trna_id = ""
    for line in maf_fh:
        if line[0] in string.ascii_letters:
            line = line.strip()
            sl = line.split()
            seq_name = sl[0]

            if not seq_name.startswith(FLOWCELL_PREFIX):
                trna_id = seq_name
                match_start = int(sl[1]) + 1
                match_end = match_start + int(sl[2])
                match_name = f"{trna_id}_{match_start}_{match_end}"
            else:
                sequence = sl[-1]
                sequence = sequence.replace("-", "")
                if match_name not in sequence_frequency_dict.keys():
                    sequence_frequency_dict[match_name] = [1, sequence]
                else:
                    sequence_frequency_dict[match_name][0] += 1


# print(sequence_frequency_dict)

for key, value in sequence_frequency_dict.items():
    matches_count = value[0]
    if matches_count >= 1:
        fasta_name = f">{key}_count{value[0]}"
        seq = value[1]
        print(fasta_name)
        print(seq)
