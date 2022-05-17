#!/usr/bin/env python 

import numpy as np
import glob
import sys
import pprint
from collections import namedtuple


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn  as sns


pp = pprint.PrettyPrinter(indent=4)


min_count = 10
min_score = 1e-03
max_match_position = 105




glob_pattern = "*.cmscan_isoall.top_hits.out"

samples_dict = {}

Isotype_data = namedtuple("Isotype_data", "coverage end_points")


def get_heatmap(numpy_file_matrix, isotypes_list, in_fn):
    title_data =  in_fn.split(".")[0]
    title = f"{title_data} heatplot of coverage (normalized by tRNA isotype mappings)"
    fig, ax = plt.subplots()
    im = ax.imshow(numpy_file_matrix)
    ax.set_yticks(np.arange(len(isotypes_list)), labels=isotypes_list)
    ax.set_xticks(np.arange(max_match_position), labels=range(1,max_match_position+1))
    ax.set_title(f"data from {title}")
    fig.tight_layout()
    plt.show()
    #plt.savefig(f"heatmap_coverage_{title_data}.png", format = "png", dpi=1200)
    #plt.close(fig) 


def seaborn_heatmap(numpy_file_matrix, isotypes_list, in_fn):
    title_data =  in_fn.split(".")[0]
    title = f"{title_data} heatplot of coverage (normalized by tRNA isotype mappings)"
    
    df = pd.DataFrame(numpy_file_matrix)
    df.index = isotypes_list
    fig, ax = plt.subplots()
    sns.heatmap(df, cmap='viridis')
    #plt.yticks = isotypes_list
    ax.set_yticks(range(0, len(isotypes_list)), labels=isotypes_list)
    ax.set_xticks(range(0, 105), labels=range(1, 105+1))
    ###plt.set_yticklabels(isotypes_list)
    plt.show()
    #plt.savefig(f"heatmap_coverage_{title_data}.png", format = "png", dpi=1200)
    #plt.close(fig) 

def normalize_isotype_data(isotype_data_input):
    last_coverage_pos = max_match_position
    
    coverage_dict   = isotype_data_input.coverage
    end_points_dict = isotype_data_input.end_points
    
    trna_isoform_coverage_vals = list(coverage_dict.values())
    total_coverage_sum = sum(trna_isoform_coverage_vals)
    for counter in range(max_match_position, 1, -1):
        if coverage_dict[counter] > 0:
            last_coverage_pos = counter
            break
    
    print(f"total_sum: {total_coverage_sum} last_coverage_pos: {last_coverage_pos}")
    print("debug_001_coverage")
    for position in range(1, last_coverage_pos+1):
        if coverage_dict[position] > 0:
            coverage_dict[position] =  100 * last_coverage_pos * (coverage_dict[position] / total_coverage_sum)
        if end_points_dict[position] > 0:
           end_points_dict[position] = 100 * last_coverage_pos * (end_points_dict[position] / total_coverage_sum)
    
    if last_coverage_pos < max_match_position:
        for  position in range(last_coverage_pos+1,max_match_position+1 ):
            end_points_dict[position] = np.nan
            coverage_dict[position] = np.nan

    print("debug coverage_dict")
    pp.pprint(coverage_dict)
    print("debug end_points_dict")
    pp.pprint(end_points_dict)
    end_points_list = list(end_points_dict.values())
    coverage_list   = list(coverage_dict.values())
    #return end_points_list
    return coverage_list

def parse_single_file(in_fn, symbol):
    print("QQQ", in_fn)
    isoforms_dict = {}
    with open(in_fn) as f:
        for line in f:
            line = line.strip()
            sl = line.split()
            isotype_cmscan = sl[0]
            fragment_count = int(sl[2].split("count")[-1])
            match_start = int(sl[5])
            match_end   = int(sl[6])
            score = float(sl[-3])
            if (fragment_count >= min_count) and (score <= min_score):
                if isotype_cmscan not in isoforms_dict.keys():
                    isoforms_dict[isotype_cmscan] = Isotype_data({}, {})
                    for i in range(1, max_match_position+1): 
                        isoforms_dict[isotype_cmscan].coverage[i] = 0
                        isoforms_dict[isotype_cmscan].end_points[i] = 0
                
                isoforms_dict[isotype_cmscan].end_points[match_start] += fragment_count
                isoforms_dict[isotype_cmscan].end_points[match_end]   += fragment_count

                for position in range(match_start, match_end+1):
                    isoforms_dict[isotype_cmscan].coverage[position] += fragment_count
                   
                

                #isoforms_dict[isotype_cmscan][match_start] += fragment_count
                #isoforms_dict[isotype_cmscan][match_end] += fragment_count
                #print(isotype_cmscan, fragment_count, match_start, match_end, score)
                #print("XXX", match_start)
        #pp.pprint(isoforms_dict)
        isotypes_one_file_list = list(isoforms_dict.keys())
        num_of_one_file_isotypes = len(isotypes_one_file_list) 
        
        numpy_file_matrix = np.arange(num_of_one_file_isotypes * max_match_position,  dtype=float)
        numpy_file_matrix = numpy_file_matrix.reshape(num_of_one_file_isotypes, max_match_position)
        numpy_file_matrix = np.zeros_like(numpy_file_matrix)
        print("foo numpy")
        print(numpy_file_matrix.shape)
        pp.pprint(numpy_file_matrix)
        for index, isotype in enumerate(isotypes_one_file_list):
            print(index, isotype)
            numpy_file_matrix[index,:] = normalize_isotype_data(isoforms_dict[isotype])
        print("normalized p1 ends")
        #pp.pprint(numpy_file_matrix)
        #get_heatmap(numpy_file_matrix, isotypes_one_file_list, in_fn)
        seaborn_heatmap(numpy_file_matrix, isotypes_one_file_list, in_fn)
    #print(in_fn)
    #print(isoforms_dict)
    return isoforms_dict

for in_fn in glob.glob(glob_pattern):
    sample_name = in_fn.split(".")[0]
    split_sample_name = sample_name.split("_")

    print(split_sample_name,sample_name, in_fn)
    if len(split_sample_name) == 3:
        treated_sample = True
        symbol = "BH4"
    else:
        treated_sample = False
        symbol = "UNTR"
    sample_grp =  f"{split_sample_name[0]}_{split_sample_name[-1]}"
    if sample_grp not in samples_dict.keys():
       samples_dict[sample_grp] = {"UNTR": {}, "BH4": {}}

    file_isoforms_dict = parse_single_file(in_fn, symbol)
    #print("XXX", sample_grp, symbol, file_isoforms_dict)
    samples_dict[sample_grp][symbol] = file_isoforms_dict


#pp.pprint(samples_dict)
#print(samples_dict)


"""
## FIXME
for sample_grp in samples_dict.keys():
    for symbol in ["UNTR", "BH4"]: 
        for trna_isoform in samples_dict[sample_grp][symbol].keys():
            trna_isoform_vals = list(samples_dict[sample_grp][symbol][trna_isoform].values())
            trna_isoform_vals = [x for x in trna_isoform_vals if x > 0 ]
            non_zeros_len = len(trna_isoform_vals)
            non_zeros_sum = sum(trna_isoform_vals)
            factor = non_zeros_sum/(non_zeros_len)

            print(sample_grp, symbol, trna_isoform)
            print(factor, non_zeros_len, non_zeros_sum)
            print(trna_isoform_vals)
            foo_list = []
            for item in trna_isoform_vals:
                foo_list.append(1000 * item/factor)
            print(foo_list)

"""

""" 
    samples_dict[sample] = []
    trna_seq_coverage = {}
    for i in range(1,101):
        trna_seq_coverage[i] = 0
    samples_dict[sample].append(trna_seq_coverage)

"""


