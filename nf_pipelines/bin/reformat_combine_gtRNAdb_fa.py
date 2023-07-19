#!/usr/bin/env pypy3

"""
combines confident and mature tRNA fasta files from gtRNAdb
adds CCA tail to all tRNAs
adds .intron (confident set) or .spld (mature) to spliced tRNA names

"""

import textwrap

from pyfaidx import Fasta

spliced_trnas = "Arg-TCT Ile-TAT Leu-CAA Tyr-GTA Tyr-ATA".split()


def split_by_len(txt: str, len: int, sep: str or None = "\n") -> str or list:
    """
    txt: str text
    len: split length (symbols per split)
    sep: separate string or None for list of strs
    """
    spl_list = [txt[i * len : i * len + len] for i in range(len(txt) // l + 1)]
    return spl_list if sep is None else sep.join(spl_list)


def add_cca_fasta(fasta_in, spliced_flag):
    def fix_seq(record, new_name):
        print(f">{new_name}")
        seq = f"{record}"
        seq = seq + "CCA"
        print("\n".join(textwrap.wrap(seq, 80)))

    in_fa = Fasta(fasta_in, as_raw=True, sequence_always_upper=True)

    seq_names_dict = {}

    for record in in_fa:
        tmp_list = f"{record.name}".split("_tRNA-")[1].split("-")
        isotype = "-".join(tmp_list[:2])
        # print(tmp_list, isotype)
        new_name = f"{record.name}".replace("Homo_sapiens_", "")
        if spliced_flag == "intron":
            if isotype in spliced_trnas:
                new_name = new_name + ".intron"
            fix_seq(record, new_name)
        elif spliced_flag == "spld":
            if isotype in spliced_trnas:
                new_name = new_name + ".spld"
                fix_seq(record, new_name)
            else:
                # new_name = new_name + ".discard"
                pass
        else:
            print("wrong flag", spliced_flag)


if __name__ == "__main__":
    mature_gtrnadb_fn = "hg38-mature-tRNAs.fa"
    confident_gtrnadb_fa = "hg38-tRNAs.fa"

    # gtrnadb_fa = sys.argv[1]
    add_cca_fasta(confident_gtrnadb_fa, "intron")
    add_cca_fasta(mature_gtrnadb_fn, "spld")
