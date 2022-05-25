#!/usr/bin/env python3

"""
runs clumpify from BBmap to reduce the size
"""

import glob
import os
import sys
import textwrap


def slurm_setup(job_name):
    """prints the SLURM shell file header"""

    shell_file_header = f"""
    #!/bin/bash"

    #SBATCH --nodes=1 
    #SBATCH --time={JOB_TIME} 
    #SBATCH --cpus-per-task={NUM_OF_THREADS}
    #SBATCH --partition={SLURM_PARTITION}

    """
    shell_file_header = textwrap.dedent(shell_file_header)

    env_setup = """
    module load java/11
    # pigz install dir
    export PATH=$HOME/soft/bin:$PATH
    """

    env_setup = textwrap.dedent(env_setup)

    print(shell_file_header)
    print(f"#SBATCH --job-name={job_name}")
    print(env_setup)


def one_file_shell(input_fastq_fn):
    """create shell file for sbatch/SLURM to run fastp"""
    # save the basename  and ending
    fastq_fn_base = input_fastq_fn.split("/")[-1]
    fastq_fn_base = fastq_fn_base.replace(".fastq.gz", "")

    out_fq_fn = f"{TOP_OUT_DIR}/{fastq_fn_base}.clump_opt_dedup.fq.gz"

    job_name = f"clump_{fastq_fn_base}"
    shell_fn = f"{job_name}.sh"

    saveout = sys.stdout

    with open(shell_fn, "w", encoding="utf-8") as output_fh:
        sys.stdout = output_fh
        slurm_setup(job_name)

        print("""echo "starting clumpify" """)
        print("date")

        command_1 = f"""
        {CLUMPIFY_EXE} \\
        dedupe=t \\
        optical=t \\
        reorder=a \\
        shortname=shrink \\
        blocksize=2048 \\
        ziplevel=8 \\
        pigz={NUM_PIGZ_THREADS} \\
        unpigz={NUM_PIGZ_THREADS} \\
        lowcomplexity=t \\
        in={input_fastq_fn} \\
        out={out_fq_fn}
        """

        command_1 = textwrap.dedent(command_1)  # .replace(8 * " ", "")

        print(command_1)

        print("\n\n")
        print("""echo "finished clumpify" """)
        print("date")

        sys.stdout = saveout
        output_fh.close()
        os.system("chmod +x {shell_fn}")


if __name__ == "__main__":

    # setup
    JOB_TIME = "05:00:00"
    SLURM_PARTITION = "mem"
    NUM_OF_THREADS = 16
    NUM_PIGZ_THREADS = 6

    ## setup clump
    INPUT_DIR = "./ORIG/"
    INPUT_PATTERN = "*.fastq.gz"

    TOP_OUT_DIR = "./"
    CLUMP_OUT_DIR = "CLUMP_20220524"
    output_dir_top = f"{TOP_OUT_DIR}/CLUMP_OUT_DIR"

    CLUMPIFY_EXE = "/scratch/dkedra/soft/progs/bbmap_current/clumpify.sh"

    for fastq_fn in glob.glob(f"{INPUT_DIR}/{INPUT_PATTERN}"):
        one_file_shell(fastq_fn)
