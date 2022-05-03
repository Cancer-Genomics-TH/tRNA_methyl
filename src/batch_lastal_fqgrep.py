#!/usr/bin/env python3

"""
mapping fqgrep-masked fastq files to GRCh38_plus_tRNAs
using lastal aligner

dir: /scratch/dkedra/proj/trna_20220426/mapping
"""

import glob
import os
import sys

input_dir = "./fq_inputs"
out_dir = "./map_results/fqgrep_masked"

# fastq masked input fastq
input_fn_pattern_1 = "*.clump_opt_dedup.fastp.fqgrep_mask.fq.gz"

patterns = [
    input_fn_pattern_1,
]
last_db = "genome_last/hg38_tRNAs/GRCh38_plus_tRNAs.last"

num_of_threads = 48

shell_file_header = f"""#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task={num_of_threads}
#SBATCH --partition=mem
"""

env_setup = """
module load gcc/11.2.0 

export PATH=/scratch/dkedra/soft/bin:$PATH

"""


for pattern in patterns:
    for fastq_fn in glob.glob(f"{input_dir}/{pattern}"):
        # save the basename  and ending
        fastq_fn_base = fastq_fn.split("/")[-1]
        fastq_fn_base = fastq_fn_base.replace(".fq.gz", "")
        out_maf_fn = f"{out_dir}/{fastq_fn_base}.lastal.hg38-tRNAs.maf.gz"

        job_name = f"lastal_{fastq_fn_base}"
        shell_fn = f"{job_name}.sh"

        saveout = sys.stdout
        output_fh = open(shell_fn, "w")
        sys.stdout = output_fh

        print(shell_file_header)
        print(f"#SBATCH --job-name={job_name}")
        print(env_setup)

        command_1 = f"""lastal -v -P{num_of_threads} -Qkeep -C2 {last_db} {fastq_fn} | last-split | gzip > {out_maf_fn} 
        """

        command_1 = command_1.replace(8 * " ", "")

        print(command_1)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system("chmod +x %s" % (shell_fn))
