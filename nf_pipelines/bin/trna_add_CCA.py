#!/usr/bin/env pypy3

import sys
import textwrap

from pyfaidx import Fasta


def split_by_len(txt: str, len: int, sep: str or None = "\n") -> str or list:
    """
    txt: str text
    len: split length (symbols per split)
    sep: separate string or None for list of strs
    """
    spl_list = [txt[i * len : i * len + len] for i in range(len(txt) // l + 1)]
    return spl_list if sep is None else sep.join(spl_list)


def add_CCA_fasta(gtrnadb_fa):
    in_fa = Fasta(gtrnadb_fa)

    seq_names_dict = {}

    for record in in_fa:
        print(f">{record.long_name}")
        seq = f"{record}"
        seq = seq + "CCA"
        print("\n".join(textwrap.wrap(seq, 80)))


if __name__ == "__main__":
    gtrnadb_fa = sys.argv[1]
    add_CCA_fasta(gtrnadb_fa)
