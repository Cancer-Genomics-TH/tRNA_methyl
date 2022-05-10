#!/usr/bin/env python3

"""
create SLURM shell scripts to
extract lines containg tRNA contigs from a gzipped .maf files

count	seq_name	match_start	match_end	seq_len	match_seq


"""

import glob
import os
import sys

input_dir = "./"
out_dir = "./"




trna_names_4rg = "/scratch/dkedra/proj/trna_20220426/mapping/trna_names_short.with_s.txt"

input_fn_pattern_1 = "*.maf.gz"

patterns = [input_fn_pattern_1, ]

num_of_threads = 4

shell_file_header = f"""#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task={num_of_threads}
#SBATCH --partition=express
"""

env_setup = """
module load gcc/11.2.0 

export PATH=/scratch/dkedra/soft/bin:$PATH

"""


for pattern in patterns:
    for maf_fn in glob.glob(f"{input_dir}/{pattern}"):
        # save the basename  and ending
        maf_fn_base = maf_fn.split("/")[-1]
        maf_fn_base = maf_fn_base.replace(".maf.gz", "")
        out_trna_tsv_fn = f"{out_dir}/{maf_fn_base}.trna.frag_cnt.tsv.gz"

        job_name = f"rg_trna_{maf_fn_base}"
        shell_fn = f"{job_name}.sh"

        saveout = sys.stdout
        output_fh = open(shell_fn, "w")
        sys.stdout = output_fh

        print(shell_file_header)
        print(f"#SBATCH --job-name={job_name}")
        print(env_setup)

        # extracts aligned lines for the tRNA contigs 
        # prints tRNA_name, match_start(+1 because MAF is zero based), match_end, 
        # adding header
        command_1 = f"cp header.txt.gz {out_trna_tsv_fn}"

        command_2 = f"""rg -z --threads {num_of_threads} -f {trna_names_4rg} {maf_fn} | awk '{{print $2, $3+1, $3+$4, $6, $7}}' | sort | uniq -c |  tr --squeeze-repeats " " | sed 's/ /\\t/g' | sed -e 's/^\\t//g' | gzip >> {out_trna_tsv_fn} 
        """

        command_1 = command_1.replace(8 * " ", "")
        command_2 = command_2.replace(8 * " ", "") 
        
        print(command_1)
        print()
        print(command_2)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system("chmod +x %s" % (shell_fn))
