#!/usr/bin/env python3
"""
* extracts the top hits from the cmscan tab file

"""


import sys



def parse_cmscan_tab(in_fn):


    with open(in_fn) as fh:
        old_query_id = ""
        for line in fh:
            if line[0] == "#":
                pass
            else: 
                #new_line = "\t".join(line.split())
                sl = line.split()
                #line = line.strip()
                try:
                    query_id = sl[2]
                except:
                    print("xxx_debug1", line)
                if query_id != old_query_id:
                    out_str = "\t".join(sl)
                    print(out_str)
                old_query_id = query_id
            

if __name__ == "__main__":
    in_fn = sys.argv[1]
    parse_cmscan_tab(in_fn)
