#!/usr/bin/env python3

"""
mapping fqgrep-masked fastq files to GRCh38_plus_tRNAs
using lastal aligner

dir: /scratch/dkedra/proj/trna_20220426/mapping
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
    module load gcc/11.2.0
    export PATH=/scratch/dkedra/soft/bin:$PATH
    """

    env_setup = textwrap.dedent(env_setup)

    print(shell_file_header)
    print(f"#SBATCH --job-name={job_name}")
    print(env_setup)


def one_file_shell(input_fastq_fn):
    """FIXME"""
    # save the basename  and ending
    fastq_fn_base = input_fastq_fn.split("/")[-1]
    fastq_fn_base = fastq_fn_base.replace(".fq.gz", "")
    out_maf_fn = f"{OUT_DIR}/{fastq_fn_base}.lastal.hg38-tRNAs.maf.gz"

    job_name = f"lastal_{fastq_fn_base}"
    shell_fn = f"{job_name}.sh"

    saveout = sys.stdout

    with open(shell_fn, "w", encoding="utf-8") as output_fh:
        sys.stdout = output_fh
        slurm_setup(job_name)

        print("""echo "starting lastal mapping" """)
        print("date")

        command_1 = f"""
        lastal -v -P{NUM_OF_THREADS} -Qkeep -C2 {LASTDB_FN} {input_fastq_fn} | last-split | gzip > {out_maf_fn}
        """

        command_1 = textwrap.dedent(command_1)

        print(command_1)
        print("\n\n")

        sys.stdout = saveout
        output_fh.close()
        os.system(f"chmod +x {shell_fn}")


if __name__ == "__main__":

    # setup
    NUM_OF_THREADS = 48
    JOB_TIME = "05:00:00"
    SLURM_PARTITION = "mem"

    INPUT_DIR = "./fq_inputs"
    OUT_DIR = "./map_results/fqgrep_masked"

    LASTDB_FN = "genome_last/hg38_tRNAs/GRCh38_plus_tRNAs.last"
    # fastq masked input fastq
    INPUT_PATTERN = "*.clump_opt_dedup.fastp.fqgrep_mask.fq.gz"

    for fastq_fn in glob.glob(f"{INPUT_DIR}/{INPUT_PATTERN}"):
        one_file_shell(fastq_fn)
