#!/usr/bin/env python3.10

import os
import sys
import textwrap

from pyfaidx import Fasta

from trnascan_out_parser import parse_out


def split_by_len(txt: str, len: int, sep: str or None = "\n") -> str or list:
    """
    txt: str text
    len: split length (symbols per split)
    sep: separate string or None for list of strs
    """
    spl_list = [txt[i * len : i * len + len] for i in range(len(txt) // l + 1)]
    return spl_list if sep is None else sep.join(spl_list)


def process_trnascan_fasta(combo_fa):
    in_fa = Fasta(combo_fa)

    seq_names_dict = {}

    for record in in_fa:
        name = record.long_name
        """
        if name.endswith("pseudogene"):
            pseudo = "pseudo"
        else:
            pseudo = ""
        """
        sl = name.split()
        trna_id = sl[0]
        if trna_id in tags_dict.keys():
            tag = tags_dict[trna_id]
            new_name = f"tRNA-{sl[3]}-{sl[4][1:-1]}.trnascan.{sl[0]}.{tag}"
        else:
            new_name = f"tRNA-{sl[3]}-{sl[4][1:-1]}.trnascan.{sl[0]}"
        ##if new_name.endswith("_"):
        ##    new_name = new_name[:-1]
        # print(new_name)
        # print(name)
        # print(sl[0], sl[3], sl[4][1:-1], pseudo)
        if new_name not in seq_names_dict.keys():
            seq_names_dict[new_name] = sl[0]
        else:
            print("ERROR", name, new_name)

    for new_name in sorted(seq_names_dict.keys()):  # .sorted:
        old_name = seq_names_dict[new_name]
        seq = f"{in_fa[old_name]}"
        print(f">{new_name}")
        seq = seq + "CCA"
        # nseq = split_by_len(seq, 40, "\n")
        print("\n".join(textwrap.wrap(seq, 80)))
        ##return seq_names_dict


def filter_chr_trnas(raw_chr_trnas_fn):
    chr_trnas = Fasta(raw_chr_trnas_fn, filt_function=lambda x: x[:2] != "MT")

    saveout = sys.stdout
    with open("filtered_chr_trnas_fa", "w") as output_fh:
        sys.stdout = output_fh
        for record in chr_trnas:
            print(f">{record.long_name}")
            print(record)

        sys.stdout = saveout
        output_fh.close()

    # mt_trnas =  Fasta(mt_trnas_fn)


if __name__ == "__main__":
    raw_chr_trnas_fn = sys.argv[1]
    mt_trnas_fn = sys.argv[2]
    chr_trnas_out = sys.argv[3]
    tags_dict = parse_out(chr_trnas_out)

    filter_chr_trnas(raw_chr_trnas_fn)
    # combine filtered chr and mt trnas
    command_1 = f"cat filtered_chr_trnas_fa {mt_trnas_fn} > raw_combined_trnas.fa"
    # print(command_1)
    os.system(command_1)
    process_trnascan_fasta("raw_combined_trnas.fa")

    # process_trnascan_fasta(chr_trnas)
