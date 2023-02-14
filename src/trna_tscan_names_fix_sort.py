#!/usr/bin/env pypy3

import sys
import textwrap
from pyfaidx import Fasta


def split_by_len(txt: str, len: int, sep: str or None='\n') -> str or list:
    """
    txt: str text
    len: split length (symbols per split)
    sep: separate string or None for list of strs
    """
    spl_list = [txt[i * len : i * len + len] for i in range(len(txt) // l + 1)]
    return spl_list if sep==None else sep.join(spl_list)



input_fn = sys.argv[1]

in_fa = Fasta(input_fn)

seq_names_dict = {}

counter = 0
for record in in_fa:
    if counter < 20:
        #print(record.long_name)
        name = record.long_name
        if name.endswith("pseudogene"):
            pseudo = "pseudo"
        else:
            pseudo = ""
        sl = name.split()
        new_name = f"tRNA-{sl[3]}-{sl[4][1:-1]}.trnascan.{sl[0]}_{pseudo}"
        if new_name.endswith("_"):
            new_name = new_name[:-1]
        print(new_name)
        print(name)
        #print(sl[0], sl[3], sl[4][1:-1], pseudo)
        if new_name not in seq_names_dict.keys():
            seq_names_dict[new_name] = sl[0]
        else:
            print("ERROR", name, new_name)
    else:
        break
    counter += 1



for new_name in sorted(seq_names_dict.keys()): #.sorted:
    old_name = seq_names_dict[new_name]
    seq = f"{in_fa[old_name]}"
    print(f">{new_name}")
    seq = seq + "aaaaaaCCA"
    #nseq = split_by_len(seq, 40, "\n")
    print("\n".join(textwrap.wrap(seq, 80)))
