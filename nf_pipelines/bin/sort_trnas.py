#!/usr/bin/env python3

import sys

from pyfaidx import Fasta

in_fa = sys.argv[1]

trnas_combo = Fasta(in_fa)

for name in sorted(trnas_combo.keys()):
    print(f">{name}")
    for line in trnas_combo[name]:
        print(line)
