#!/usr/bin/env python

"""
1. parse pre-processed cmscan from Infernal tabuleted outputs cmscan_isoall.top_hits.out
2. get the fragment positions and counts
3. normalize the counts per tRNA isoform
4. plot the coverage heat maps, one per sample

"""


import glob
import logging as log
import pickle
import pprint
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

pp = pprint.PrettyPrinter(indent=4)


samples_dict = {}

Isotype_data = namedtuple("Isotype_data", "coverage end_points")


def seaborn_heatmap(numpy_file_matrix, isotypes_list, in_fn):
    """FIXME"""
    title_data = in_fn.split(".")[0]
    title = f"Normalized Coverage Heat Plot of {title_data}"
    ## (normalized by tRNA isotype mappings)"

    df = pd.DataFrame(numpy_file_matrix)
    df.index = isotypes_list
    log.info(f"""df_hape:, {df.shape}""")

    fig, ax = plt.subplots()

    sns.heatmap(df, cmap="Reds", linewidth=1, linecolor="w", square=True)

    ax.set_yticks(
        range(0, len(isotypes_list)), labels=isotypes_list, fontname="monospace"
    )
    ax.set_xticks(range(0, 105), labels=range(1, 105 + 1))

    ax.set_xlabel("tRNA positions")
    ax.set_ylabel("isoforms")

    # heatmap.subplots_adjust(top=.8)
    # FIXME: check https://www.statology.org/seaborn-title/
    ax.set_title(f"{title}")
    # plt.text(x=4.7, y=4.7, s='Coverage Heat Plot', fontsize=16, weight='bold')
    # plt.text(x=4.7, y=4.6, s=f'{in_fn}', fontsize=8, alpha=0.75)

    # plt.set_yticklabels(isotypes_list)
    plt.show()
    # FIXME: check the resolution
    # PNGs or as below SVGs with plots are not OK
    # plt.savefig(f"heatmap_coverage_{title_data}.svg", format = "svg")
    # plt.close(fig)


def normalize_isotype_data(isotype_data_input, isotype):
    """FIXME"""

    last_coverage_pos = MAX_MATCH_POSITION

    coverage_dict = isotype_data_input.coverage
    # end_points_dict = isotype_data_input.end_points

    trna_isoform_coverage_vals = list(coverage_dict.values())
    total_coverage_sum = sum(trna_isoform_coverage_vals)
    for counter in range(MAX_MATCH_POSITION, 1, -1):
        if coverage_dict[counter] > 0:
            last_coverage_pos = counter
            break

    log.info(f"total_sum: {total_coverage_sum} last_coverage_pos: {last_coverage_pos}")

    if isotype in trna_sizes_dict.keys():
        isotype_size = trna_sizes_dict[isotype]
        log.info(
            f"isotype: {isotype} iso_size: {isotype_size} last_coverage_pos: {last_coverage_pos}"
        )
        if isotype_size > last_coverage_pos:
            last_coverage_pos = isotype_size

    for position in range(1, last_coverage_pos + 1):
        if coverage_dict[position] > 0:
            coverage_dict[position] = (
                1000 * coverage_dict[position] / total_coverage_sum
            )

    if last_coverage_pos < MAX_MATCH_POSITION:
        for position in range(last_coverage_pos + 1, MAX_MATCH_POSITION + 1):
            coverage_dict[position] = np.nan

    log.debug("debug coverage_dict")
    log.debug(coverage_dict)

    coverage_list = list(coverage_dict.values())
    # return end_points_list
    return coverage_list


def parse_single_file(in_fn, symbol):
    """FIXME"""
    isoforms_dict = {}
    with open(in_fn) as fh:
        for line in fh:
            line = line.strip()
            sl = line.split()
            isotype_cmscan = sl[0]
            fragment_count = int(sl[2].split("count")[-1])
            match_start = int(sl[5])
            match_end = int(sl[6])
            score = float(sl[-3])
            if (fragment_count >= MIN_COUNT) and (score <= MIN_SCORE):
                if isotype_cmscan not in isoforms_dict.keys():
                    isoforms_dict[isotype_cmscan] = Isotype_data({}, {})
                    for i in range(1, MAX_MATCH_POSITION + 1):
                        isoforms_dict[isotype_cmscan].coverage[i] = 0

                for position in range(match_start, match_end + 1):
                    isoforms_dict[isotype_cmscan].coverage[position] += fragment_count

        isotypes_one_file_list = list(isoforms_dict.keys())
        num_of_one_file_isotypes = len(isotypes_one_file_list)

        numpy_file_matrix = np.arange(
            num_of_one_file_isotypes * MAX_MATCH_POSITION, dtype=float
        )
        numpy_file_matrix = numpy_file_matrix.reshape(
            num_of_one_file_isotypes, MAX_MATCH_POSITION
        )
        numpy_file_matrix = np.zeros_like(numpy_file_matrix)
        log.info(f"""numpy_file_matrix_shape: {numpy_file_matrix.shape}""")
        log.debug(numpy_file_matrix)
        size_sorted_isoforms_one_file_list = []
        for sorted_tuple in isoform_sizes_sorted:
            if sorted_tuple[0] in isotypes_one_file_list:
                size_sorted_isoforms_one_file_list.append(sorted_tuple[0])

        isotypes_one_file_list = size_sorted_isoforms_one_file_list
        for index, isotype in enumerate(isotypes_one_file_list):
            log.debug(f"""filename: {in_fn} isotypes:""")
            log.debug(f"""{index}, {isotype}""")
            numpy_file_matrix[index, :] = normalize_isotype_data(
                isoforms_dict[isotype], isotype
            )
        ##print("normalized p1 ends")
        # pp.pprint(numpy_file_matrix)
        # get_heatmap(numpy_file_matrix, isotypes_one_file_list, in_fn)
        seaborn_heatmap(numpy_file_matrix, isotypes_one_file_list, in_fn)
    return isoforms_dict


def process_files(glob_pattern):
    """FIXME"""
    for in_fn in glob.glob(glob_pattern):
        sample_name = in_fn.split(".")[0]
        split_sample_name = sample_name.split("_")

        log.info(f"""{split_sample_name}, {sample_name}, {in_fn}""")
        if len(split_sample_name) == 3:
            treated_sample = True
            symbol = "BH4"
        else:
            treated_sample = False
            symbol = "UNTR"
        sample_grp = f"{split_sample_name[0]}_{split_sample_name[-1]}"
        if sample_grp not in samples_dict.keys():
            samples_dict[sample_grp] = {"UNTR": {}, "BH4": {}}

        file_isoforms_dict = parse_single_file(in_fn, symbol)
        samples_dict[sample_grp][symbol] = file_isoforms_dict


if __name__ == "__main__":
    # logging setup
    # debug_mode = False
    DEBUG_MODE = True

    ### set the logging level ###

    if DEBUG_MODE is False:
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s:%(levelname)s:%(message)s",
        )
    else:
        ## logged to file DEBUG ##
        log.basicConfig(
            filename="parse_and_plot.log",
            level=log.DEBUG,
            format="%(asctime)s:%(levelname)s:%(message)s",
        )

    # settings
    PICKLE_FN = "trna_sizes.pickle"
    MIN_COUNT = 10
    MIN_SCORE = 1e-03
    MAX_MATCH_POSITION = 105
    GLOB_PATTERN = "*.cmscan_isoall.top_hits.out"

    # get the tRNA idoforms sizes
    with open(PICKLE_FN, "rb") as f:
        trna_sizes_dict = pickle.load(f)

    isoform_sizes_sorted = sorted(trna_sizes_dict.items(), key=lambda x: x[1])
    # pp.pprint(test_sizes)
    log.debug(trna_sizes_dict)

    process_files(GLOB_PATTERN)
