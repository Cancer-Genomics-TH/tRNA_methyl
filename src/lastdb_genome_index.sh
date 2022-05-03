#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=24
#SBATCH --partition=mem

module load gcc/11.2.0

export PATH=/scratch/dkedra/soft/bin:$PATH

lastdb -P24 -uNEAR GRCh38_plus_tRNAs.last GRCh38_plus_tRNAs.fa

## dir: /scratch/dkedra/proj/trna_20220426/mapping/genome_last/hg38_tRNAs
