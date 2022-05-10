#!/home/darked89/soft/progs/pypy3.8_7.3.7/bin/pypy3

"""
uses pypy3 for speed

script to transform 2 lines of maf aligment to a fasta file
with a fasta names containing:

tRNA contig name and the number of counts
the sequence is from the NGS read match (gaps introduced during the mapping are removed )

"""


import sys

in_fn = sys.argv[1]

sequence_frequency_dict = {}

with open(in_fn) as f:
    old_trna_id = ""
    for line in f:
        if line[0] == "s":
            line = line.strip()
            sl = line[2:].split()
            seq_name = sl[0]
            #print(seq_name)
            if seq_name[4:8] == "tRNA":
                trna_id = seq_name
                match_start = int(sl[1]) + 1
                match_end = match_start + int(sl[2])
                match_name = f"{trna_id}_{match_start}_{match_end}"
            else:
                sequence = sl[-1]
                sequence = sequence.replace("-","")
                if match_name not in sequence_frequency_dict.keys():
                    sequence_frequency_dict[match_name] = [1, sequence]
                else:
                    sequence_frequency_dict[match_name][0] += 1
                #print(f"{match_name}\t{sequence}")

#print(sequence_frequency_dict)

for key, value in sequence_frequency_dict.items():
    matches_count = value[0]
    if matches_count > 9:
        fasta_name = f">{key}_count{value[0]}"
        seq = value[1]
        print(fasta_name)
        print(seq)