#!/usr/bin/env python3

"""
runs clumpify from BBmap to reduce the size
"""

import glob
import os
import sys

input_dir = "./ORIG/"
fastq_pattern = "*.fastq.gz"

output_dir_top = "./"
out_dir = f"{output_dir_top}/CLUMP_20220426"

clumpify_exe = "/scratch/dkedra/soft/progs/bbmap_current/clumpify.sh"
num_of_threads = 16
num_pigz_threads = 6


shell_file_header = f"""#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task={num_of_threads}
#SBATCH --partition=mem
"""

env_setup = """
# module load gcc/11.2.0
module load java/11

export PATH=$HOME/soft/bin:$PATH
"""


for fastq_fn in glob.glob(f"{input_dir}/{fastq_pattern}"):
    # save the basename  and ending
    fastq_fn_base = fastq_fn.split("/")[-1]
    fastq_fn_base = fastq_fn_base.replace(".fastq.gz", "")

    out_fq_fn = f"{out_dir}/{fastq_fn_base}.clump_opt_dedup.fq.gz"

    job_name = f"clump_{fastq_fn_base}"
    shell_fn = f"{job_name}.sh"

    saveout = sys.stdout
    output_fh = open(shell_fn, "w")
    sys.stdout = output_fh

    print(shell_file_header)
    print(f"#SBATCH --job-name=job_name")
    print(env_setup)

    command_1 = f"""{clumpify_exe} \\
        dedupe=t \\
        optical=t \\
        reorder=a \\
        shortname=shrink \\
        blocksize=2048 \\
        ziplevel=8 \\
        pigz={num_pigz_threads} \\
        unpigz={num_pigz_threads} \\
        lowcomplexity=t \\
        in={fastq_fn} \\
        out={out_fq_fn}
    """

    command_1 = command_1.replace(8 * " ", "")

    print(command_1)
    print("\n\n")

    sys.stdout = saveout
    output_fh.close()
    os.system("chmod +x %s" % (shell_fn))
