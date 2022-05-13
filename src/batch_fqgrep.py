#!/usr/bin/env python3

"""
filtering and quality checking fastq files using fastp 

"""

import glob
import os
import sys

out_dir = "./FQGREP_20220426"
input_dir = "./FASTP_20220426"

input_fn_pattern_1 = "*.fq.gz"
patterns = [
    input_fn_pattern_1,
]

primers_seq_list = [
    "TGGAATTCTCGGGTGCCAAGGC",
    "TGGAATTCTCGGGTGCCAAGG",
    "TGGAATTCTCGGGTGCCAAG",
    "TGGAATTCTCGGGTGCCAA",
    "TGGAATTCTCGGGTGCCA",
    "GTTCAGAGTTCTACAGTCCGACGATC",
    "GTTCAGAGTTCTACAGTCCGACGAT",
    "GTTCAGAGTTCTACAGTCCGACGA",
    "GTTCAGAGTTCTACAGTCCGACG",
    "GTTCAGAGTTCTACAGTCCGAC",
    "GCCTTGGCACCCGAGAATTCCA",
    "GATCGTCGGACTGTAGAACTCTGAAC",
]


regex_for_fqgrep = "|".join(primers_seq_list)

# debug
# print(regex_for_fqgrep)

num_of_threads = 1

shell_file_header = f"""#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task={num_of_threads}
#SBATCH --partition=express
"""

env_setup = """
module load gcc/11.2.0 

export LD_LIBRARY_PATH=/scratch/dkedra/soft/lib:$LD_LIBRARY_PATH
export PATH=/scratch/dkedra/soft/bin:$PATH

"""


for pattern in patterns:
    for fastq_fn in glob.glob(f"{input_dir}/{pattern}"):
        # save the basename  and ending
        fastq_fn_base = fastq_fn.split("/")[-1]
        fastq_fn_base = fastq_fn_base.replace(".fq.gz", "")
        out_report_fn = f"{out_dir}/{fastq_fn_base}.fqgrep_report_pat2.out"

        job_name = f"fqgrep2_{fastq_fn_base}"
        shell_fn = f"{job_name}.sh"

        saveout = sys.stdout
        output_fh = open(shell_fn, "w")
        sys.stdout = output_fh

        print(shell_file_header)
        print(f"#SBATCH --job-name=job_name")
        print(env_setup)

        command_1 = f"""fqgrep -r -a -e -p '{regex_for_fqgrep}' -o {out_report_fn}  {fastq_fn} 
        """

        command_1 = command_1.replace(8 * " ", "")

        print(command_1)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system("chmod +x %s" % (shell_fn))
