#!/usr/bin/env python3.10

import os
import sys
import textwrap

from pyfaidx import Fasta

from trnascan_out_parser import parse_out

old_tags_dict = {
    "1.trna5": "pseudo",
    "1.trna7": "IPD",
    "1.trna8": "pseudo",
    "1.trna13": "pseudo",
    "1.trna14": "pseudo",
    "1.trna15": "intron",
    "1.trna17": "pseudo",
    "1.trna23": "IPD",
    "1.trna24": "pseudo_IPD",
    "1.trna26": "pseudo",
    "1.trna28": "pseudo",
    "1.trna37": "pseudo",
    "1.trna39": "pseudo",
    "1.trna42": "pseudo",
    "1.trna46": "pseudo",
    "1.trna48": "pseudo_IPD",
    "1.trna97": "IPD",
    "1.trna100": "IPD",
    "1.trna102": "pseudo",
    "1.trna107": "pseudo",
    "1.trna108": "IPD",
    "1.trna109": "intron",
    "1.trna111": "pseudo_IPD",
    "1.trna117": "pseudo_IPD",
    "1.trna119": "pseudo_IPD",
    "1.trna194": "pseudo_IPD",
    "1.trna195": "pseudo_IPD",
    "1.trna197": "pseudo",
    "1.trna202": "pseudo",
    "1.trna204": "pseudo_IPD",
    "1.trna205": "pseudo",
    "1.trna210": "pseudo",
    "1.trna216": "pseudo_IPD",
    "1.trna219": "IPD",
    "1.trna220": "pseudo",
    "1.trna222": "pseudo",
    "1.trna224": "pseudo",
    "1.trna225": "pseudo",
    "1.trna230": "pseudo_IPD",
    "1.trna241": "pseudo",
    "1.trna244": "pseudo",
    "1.trna250": "pseudo",
    "2.trna1": "IPD",
    "2.trna2": "intron",
    "2.trna4": "IPD",
    "2.trna5": "intron",
    "2.trna6": "pseudo_IPD",
    "2.trna8": "pseudo",
    "2.trna9": "pseudo_IPD",
    "2.trna10": "pseudo",
    "2.trna12": "pseudo",
    "2.trna13": "intron",
    "2.trna14": "pseudo_IPD",
    "2.trna15": "pseudo",
    "2.trna16": "pseudo_IPD",
    "2.trna17": "pseudo",
    "2.trna19": "IPD",
    "2.trna21": "pseudo_IPD",
    "2.trna22": "pseudo_IPD",
    "2.trna24": "pseudo",
    "3.trna2": "pseudo",
    "3.trna8": "pseudo",
    "3.trna9": "pseudo_IPD",
    "3.trna10": "pseudo_trunc_start",
    "3.trna12": "pseudo",
    "4.trna2": "pseudo_IPD",
    "5.trna1": "pseudo_intron",
    "5.trna17": "pseudo_intron",
    "5.trna21": "pseudo_IPD",
    "6.trna14": "intron",
    "6.trna15": "intron",
    "6.trna16": "intron",
    "6.trna17": "intron",
    "6.trna30": "intron",
    "6.trna41": "pseudo",
    "6.trna43": "IPD",
    "6.trna48": "IPD",
    "6.trna54": "intron",
    "6.trna57": "intron",
    "6.trna58": "pseudo_intron",
    "6.trna64": "pseudo",
    "6.trna65": "intron",
    "6.trna76": "intron",
    "6.trna82": "pseudo",
    "6.trna84": "pseudo",
    "6.trna89": "pseudo_IPD",
    "6.trna91": "pseudo",
    "6.trna106": "intron",
    "6.trna109": "IPD",
    "6.trna118": "pseudo",
    "6.trna123": "pseudo",
    "6.trna125": "pseudo",
    "6.trna129": "pseudo",
    "6.trna147": "intron",
    "6.trna148": "intron",
    "6.trna158": "pseudo_IPD",
    "7.trna1": "pseudo_trunc_start",
    "7.trna2": "pseudo",
    "7.trna6": "pseudo",
    "7.trna7": "pseudo_IPD",
    "7.trna12": "IPD",
    "7.trna15": "IPD",
    "7.trna29": "pseudo_IPD",
    "7.trna31": "pseudo_trunc_start",
    "7.trna32": "pseudo",
    "8.trna3": "intron",
    "8.trna4": "intron",
    "8.trna7": "pseudo",
    "8.trna10": "IPD",
    "8.trna11": "pseudo",
    "8.trna12": "intron",
    "9.trna1": "pseudo",
    "9.trna2": "pseudo",
    "9.trna3": "pseudo",
    "9.trna5": "intron",
    "9.trna6": "pseudo",
    "9.trna9": "pseudo",
    "10.trna3": "pseudo",
    "10.trna5": "pseudo_IPD",
    "11.trna1": "IPD",
    "11.trna3": "intron",
    "11.trna8": "pseudo",
    "11.trna10": "pseudo_IPD",
    "11.trna19": "pseudo",
    "11.trna20": "pseudo_IPD",
    "12.trna9": "pseudo",
    "12.trna14": "pseudo_IPD",
    "12.trna15": "pseudo_IPD",
    "13.trna1": "pseudo_IPD",
    "14.trna5": "intron",
    "14.trna13": "pseudo",
    "14.trna14": "pseudo",
    "14.trna15": "intron",
    "14.trna16": "intron",
    "14.trna17": "intron",
    "14.trna18": "intron",
    "14.trna19": "intron",
    "15.trna5": "pseudo_trunc_start",
    "15.trna6": "pseudo",
    "16.trna5": "intron",
    "16.trna14": "pseudo",
    "16.trna15": "pseudo_IPD",
    "16.trna23": "pseudo_IPD",
    "16.trna32": "pseudo",
    "17.trna3": "intron",
    "17.trna15": "IPD",
    "17.trna20": "pseudo_IPD",
    "17.trna23": "IPD",
    "17.trna24": "pseudo",
    "17.trna29": "pseudo",
    "17.trna31": "pseudo_IPD",
    "17.trna38": "pseudo_IPD",
    "18.trna1": "pseudo",
    "19.trna3": "IPD",
    "19.trna6": "trunc_start",
    "19.trna7": "pseudo",
    "19.trna9": "pseudo",
    "19.trna10": "intron",
    "20.trna1": "pseudo",
    "20.trna2": "pseudo",
    "22.trna2": "pseudo",
    "X.trna1": "pseudo_IPD",
    "X.trna2": "pseudo_IPD",
    "X.trna3": "IPD",
    "Y.trna1": "IPD",
}


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
