#!/usr/bin/env python

import sys


in_fn = sys.argv[1]


with open(in_fn) as f:
    old_query_id = ""
    for line in f:
        line = line.strip()

        if line[0] == "#":
            #print(line)
            pass
        else:
            sl = line.split()
            try:
                query_id = sl[2]
            except:
                print("xxx_debug1", line)
            if query_id != old_query_id:
                print(line)
            old_query_id = query_id

