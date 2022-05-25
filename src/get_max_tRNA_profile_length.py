#!/usr/bin/env python3

"""
get the right most position from the cmscan tab results
dump dictionary of tRNA isoforms names : lengths as a pickle
this pickle/dictionary is needed by: parse_and_plot_heatmap.py
"""

import pickle
import pprint
import sys

pp = pprint.PrettyPrinter(indent=4)


cmscan_tabout_fn = sys.argv[1]

sizes_dict = {}

with open(cmscan_tabout_fn, encoding="utf-8") as tabout_fh:
    for line in tabout_fh:
        line = line.strip()
        sl = line.split()
        trna_name, match_start, match_end = sl
        match_end = int(match_end)
        if trna_name not in sizes_dict.keys():
            sizes_dict[trna_name] = match_end
        else:
            if match_end > sizes_dict[trna_name]:
                sizes_dict[trna_name] = match_end


# debug: pp.pprint(sizes_dict)
with open("trna_sizes.pickle", "wb") as f:
    pickle.dump(sizes_dict, f, pickle.HIGHEST_PROTOCOL)
