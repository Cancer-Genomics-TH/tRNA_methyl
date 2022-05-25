#!/usr/bin/env pypy3

"""
creates a fastq file out of a fqgrep report

input: gzipped fqgrep report
output: fastq fn printed on the STDOUT

usage:

parse_fqgrep_report.py some_report.gz > output.fastq

"""


import gzip
import sys


def parse_list(fqgrep_columns):
    """transform a list from one report row into fastq for one read"""
    read_name = fqgrep_columns[0]
    read_seq = fqgrep_columns[-2]
    read_qual = fqgrep_columns[-1]
    match = fqgrep_columns[-3]

    if match == "*":
        out_seq = read_seq
        out_qual = read_qual
        # debug: print(f"no_match: {fqgrep_columns}")
    else:
        start = int(fqgrep_columns[-5])
        end = int(fqgrep_columns[-4])
        # print(start, end)
        seq_mask = "N" * (end - start)
        qual_mask = "B" * (end - start)
        out_seq = f"{read_seq[:start]}{seq_mask}{read_seq[end:]}"
        out_qual = f"{read_qual[:start]}{qual_mask}{read_qual[end:]}"
        # debug: print(out_seq)
        # debug: print(f"got_match: {fqgrep_columns}")
    print(f"@{read_name}")
    print(out_seq)
    print("+")
    print(out_qual)


def report_to_fastq(input_fn):
    """parse the gzipped fqgrep report to masked fastq"""
    with gzip.open(input_fn, mode="rt") as fh:
        counter = 0
        for line in fh:
            if counter == 0:
                pass
            else:
                line = line.strip()
                sl = line.split()
                parse_list(sl)
            counter += 1


if __name__ == "__main__":
    input_report_fn = sys.argv[1]
    report_to_fastq(input_report_fn)
