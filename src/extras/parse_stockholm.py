#!/usr/bin/env pypy3.9

"""
find the consensus position of an anticodon from a Stockholm tRNA alignment     

outputs fasta file:

>tRNA-Cys-GCA_prof_anticodon
GGGNNNNTAGNTNANNNGGNNAGAGCANTTGACTGCA
>tRNA-Glu-CTC_prof_anticodon
TCCCTGGTGGTCTAGTGGTNAGGATTCGGCGCTC
>tRNA-Glu-TTC_prof_anticodon
TCCCNNNTGGTCTAGNGGNNAGGATTCNNNGNTTTC

tag: auxiliary script
"""


import logging as log
import re
import sys

from Bio import AlignIO

# set low, i.e 5 for testing/debug
MAX_NUM_ALIGN = 100

if len(sys.argv) > 1:
    stockholm_fn = sys.argv[1]
else:
    stockholm_fn = "gtrna_align.isotypes.all.sto"

with open(stockholm_fn) as myFile:
    counter = 0
    isoforms_dict = {}
    for alignment in AlignIO.StockholmIO.StockholmIterator(myFile):
        if counter < MAX_NUM_ALIGN:
            ali_lenght = alignment.get_alignment_length()
            log.debug(f"ali_lenght: {ali_lenght}")

            ali_num_records = len(alignment._records)
            log.debug(ali_num_records)

            isoform_raw_name = alignment[0].id
            isoform_parts_list = isoform_raw_name.split("-")[:3]
            isoform_name = "-".join(isoform_parts_list)
            isoform_anticodon = isoform_parts_list[-1]
            log.debug(isoform_name)
            log.debug(isoform_anticodon)
            log.debug(isoform_raw_name)

            frequency_dict = {}
            for record in alignment:
                for pos, base in enumerate(record.seq):
                    if pos not in frequency_dict.keys():
                        frequency_dict[pos] = {}
                    if base not in frequency_dict[pos].keys():
                        frequency_dict[pos][base] = 1
                    else:
                        frequency_dict[pos][base] += 1

                log.debug(record.id)
                log.debug(record.seq)

            annot_list = []
            for pos in frequency_dict.keys():
                bases = list(frequency_dict[pos])
                if len(bases) == 1:
                    if bases[0] == "-":
                        pass
                    else:
                        annot_list.append(bases[0])
                else:
                    annot_list.append("N")

            log.debug(f"profile_len: {len(annot_list)}")
            log.debug("annot_list: {annot_list}")

            isoform_consensus_seq = ""
            for base in annot_list:
                if base in "ACTG":
                    isoform_consensus_seq = f"{isoform_consensus_seq}{base}"
                else:
                    isoform_consensus_seq = f"{isoform_consensus_seq}N"

            anticodons_list = [
                m.start() for m in re.finditer(isoform_anticodon, isoform_consensus_seq)
            ]

            anti_candidates = [x for x in anticodons_list if x in range(31, 40)]
            anti_position = anti_candidates[0]
            print(f">{isoform_name}_prof_anticodon")
            print(isoform_consensus_seq[: anti_position + 3])
            # print(f"{isoform_name}\t{isoform_anticodon}\t{anti_candidates}\t{len(isoform_consensus_seq)}\t{isoform_consensus_seq}")

        counter += 1
        tmp_list_len = len(annot_list)
        if tmp_list_len < 105:
            for i in range(0, 105 - tmp_list_len):
                annot_list.append("")
        # print(annot_list)
        # print(len(annot_list))