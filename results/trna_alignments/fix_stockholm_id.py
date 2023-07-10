#!/usr/bin/env python3

import glob

pattern = "*.out/results/result.stk"

file_list = sorted(glob.glob(pattern))

#print(file_list)

for fn in file_list:
    isoform = fn.split(".out")[0]
    #print(fn, isoform)
    with open(fn) as fh:
        text = fh.read()
        lines = text.split("\n")
        for header_line in lines[:3]:
            print(header_line)
        out_str = f"#=GF ID {isoform}"
        print(out_str)
        for alig_line in lines[3:]:
            print(alig_line)
