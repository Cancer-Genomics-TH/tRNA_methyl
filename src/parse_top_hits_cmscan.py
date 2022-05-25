#!/usr/bin/env python

"""
keep just the top hit from the tabular output
of the cmscan from Infernal

usage:

parse_top_hits_cmscan.py  cmscan.tabular.out > cmscan.tabular.top_hits
"""

import sys

in_fn = sys.argv[1]

with open(in_fn, encoding="utf-8") as in_fh:
    old_query_id = ""
    for line in in_fh:
        line = line.strip()
        if line[0] == "#":
            pass
        else:
            sl = line.split()
            try:
                query_id = sl[2]
            except IndexError:
                print("IndexError", line)

            if query_id != old_query_id:
                print(line)
            old_query_id = query_id
