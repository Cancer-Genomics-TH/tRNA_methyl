#!/usr/bin/env python3

"""
create SLURM shell scripts to
extract lines containg tRNA contigs from a gzipped .maf files

FIXME: header.txt.gz
FIXME: "/scratch/dkedra/proj/trna_20220426/mapping/trna_names_short.with_s.txt"

count	seq_name	match_start	match_end	seq_len	match_seq

"""

import glob
import os
import sys
import textwrap


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

    export PATH=/scratch/dkedra/soft/bin:$PATH
    """

    env_setup = textwrap.dedent(env_setup)

    print(shell_file_header)
    print(f"#SBATCH --job-name={job_name}")
    print(env_setup)


def one_file_shell(in_maf_fn):
    """creates shell for SLURM for parsing one maf file"""
    # save the basename  and ending
    # save the basename  and ending
    maf_fn_base = in_maf_fn.split("/")[-1]
    maf_fn_base = maf_fn_base.replace(".maf.gz", "")
    out_trna_tsv_fn = f"{OUT_DIR}/{maf_fn_base}.trna.frag_cnt.tsv.gz"

    job_name = f"rg_trna_{maf_fn_base}"
    shell_fn = f"{job_name}.sh"

    saveout = sys.stdout
    with open(shell_fn, "w", encoding="utf-8") as output_fh:
        sys.stdout = output_fh
        slurm_setup(job_name)

        # extracts aligned lines for the tRNA contigs
        # prints tRNA_name, match_start(+1 because MAF is zero based), match_end,
        # adding header
        command_1 = f"""
        cp header.txt.gz {out_trna_tsv_fn}
        """

        command_2 = f"""
        rg -z --threads {NUM_OF_THREADS} -f {TRNA_NAMES} {in_maf_fn}
        """
        command_2 = f"""
        {command_2} | \\
        awk '{{print $2, $3+1, $3+$4, $6, $7}}' 
        """

        command_2 = f"""
        {command_2} | \\
        sort | uniq -c |  tr --squeeze-repeats ' ' 
        """

        command_2 = f"""
        {command_2} | \\
        sed 's/ /\\t/g' | sed -e 's/^\\t//g' | gzip >> {out_trna_tsv_fn} 
        """

        command_1 = textwrap.dedent(command_1)
        command_2 = textwrap.dedent(command_2)

        print(command_1)
        print()
        print(command_2)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system("chmod +x {shell_fn}")


if __name__ == "__main__":

    INPUT_DIR = "./"
    OUT_DIR = "./"
    SLURM_PARTITION = "express"
    JOB_TIME = "00:30:00"
    NUM_OF_THREADS = 4
    # trna names for ripgrep with starting "s trna_name"
    TRNA_NAMES = (
        "/scratch/dkedra/proj/trna_20220426/mapping/trna_names_short.with_s.txt"
    )

    INPUT_PATTERN = "*.maf.gz"

    for maf_fn in glob.glob(f"{INPUT_DIR}/{INPUT_PATTERN}"):
        one_file_shell(maf_fn)
