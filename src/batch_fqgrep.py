#!/usr/bin/env python3

"""
creating shell scripts for sbatch/SLURM
using fqgrep to detect seq kit primers in the fastq files

input: fastq as fq.gz
output: fqgrep report
"""

import glob
import os
import sys
import textwrap


def join_primers():
    """creates a regex string for fqgrep"""
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

    my_regex = "|".join(primers_seq_list)
    return my_regex


def slurm_setup(job_name):
    """prints the SLURM shell file header"""

    print("/usr/bin/bash") 

    shell_file_header = f"""

    #SBATCH --nodes=1 
    #SBATCH --time={JOB_TIME} 
    #SBATCH --cpus-per-task={NUM_OF_THREADS}
    #SBATCH --partition={SLURM_PARTITION}

    """
    shell_file_header = textwrap.dedent(shell_file_header)

    env_setup = """
    module load gcc/11.2.0
    
    export LD_LIBRARY_PATH=/scratch/dkedra/soft/lib:$LD_LIBRARY_PATH
    export PATH=/scratch/dkedra/soft/bin:$PATH
    """

    env_setup = textwrap.dedent(env_setup)

    print(shell_file_header)
    print(f"#SBATCH --job-name={job_name}")
    print(env_setup)


def one_file_shell(input_fastq_fn):
    """create shell for sbatch to run fqgrep on a single fastq input"""
    # save the basename  and ending
    fastq_fn_base = input_fastq_fn.split("/")[-1]
    fastq_fn_base = fastq_fn_base.replace(".fq.gz", "")
    out_report_fn = f"{OUT_DIR}/{fastq_fn_base}.fqgrep_report_pat2.out"

    job_name = f"fqgrep2_{fastq_fn_base}"
    shell_fn = f"{job_name}.sh"

    saveout = sys.stdout

    with open(shell_fn, "w", encoding="utf-8") as output_fh:
        sys.stdout = output_fh
        slurm_setup(job_name)

        print("""\necho "starting fqgrep" """)
        command_1 = f"""
        fqgrep -r -a -e -p '{regex_for_fqgrep}' -o {out_report_fn}  {input_fastq_fn} 
        """
        command_1 = textwrap.dedent(command_1)

        print(command_1)
        print("""\necho "finished fqgrep" """)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system(f"chmod +x {shell_fn}")


if __name__ == "__main__":
    INPUT_DIR = "./FASTP_20220426"
    OUT_DIR = "./FQGREP_20220426"

    SLURM_PARTITION = "express"
    NUM_OF_THREADS = 1
    JOB_TIME = "00:50:00"

    regex_for_fqgrep = join_primers()

    INPUT_FASTQ_PATTERN = "*.fq.gz"
    for fastq_fn in glob.glob(f"{INPUT_DIR}/{INPUT_FASTQ_PATTERN}"):
        one_file_shell(fastq_fn)
