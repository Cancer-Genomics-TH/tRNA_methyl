#!/usr/bin/env pypy3


"""
env:
/scratch/dkedra/.conda/envs/python_310/

* runs mlocarna on spliced tRNA fasta files
* creates sto file for single sequence:  tRNA_Tyr-ATA.fa


"""

import glob
import sys

import executor

exception_isotype_fn = "tRNA_Tyr-ATA.fa"


tyr_ata_sto = """ 
# STOCKHOLM 1.0
#=GF CC hand made
#=GF ID spld_Tyr-ATA

spld-Tyr-ATA-1-1   CCTTCAATAGTTCAGCTGGTAGAGCAGAGGACTATAGGTCCTTAGGTTGCTGGTTCGATTCCAGCTTGAAGGA
//

"""


def run_locarna(multi_fasta_fn):
    
    command_1 = f"mlocarna  -q --stockholm --threads 6 -LP {multi_fasta_fn}"
    executor.execute(command_1)

def fix_stockholm_locarna(sto_fn):
    isoform = sto_fn.split(".out")[0].split("/")[2]
    #debug
    #print(isoform)
    with open(sto_fn) as fh:
        text = fh.read()
        lines = text.split("\n")
        for header_line in lines[:3]:
            print(header_line)
        out_str = f'''#=GF ID {isoform.replace("tRNA", "spld")}'''
        print(out_str)
        for alig_line in lines[3:]:
            print(alig_line)




""" 

pattern = "*.out/results/result.stk"

file_list = sorted(glob.glob(pattern))

#print(file_list)

for fn in file_list:
    isoform = fn.split(".out")[0]
    #print(fn, isoform)
    with open(fn) as fh:
        text = fh.read()
        lines = text.split("\n")
        for header_line in lines[:3]:
            print(header_line)
        out_str = f"#=GF ID {isoform}"
        print(out_str)
        for alig_line in lines[3:]:
            print(alig_line)
"""



if __name__ == "__main__":
    spliced_fasta_dir = sys.argv[1]

    tyr_ata_fn = f"./{spliced_fasta_dir}/{exception_isotype_fn}"
    
    fasta_fn_list = glob.glob(f"./{spliced_fasta_dir}/tRNA_*.fa")
   
    index_for_del = fasta_fn_list.index(tyr_ata_fn)
    del fasta_fn_list[index_for_del]
    for fasta_fn in fasta_fn_list:
        run_locarna(fasta_fn)
    
    locarna_sto_list = glob.glob(f"./{spliced_fasta_dir}/tRNA_*.out/results/result.stk")
    for raw_sto_fn in locarna_sto_list:
        fix_stockholm_locarna(raw_sto_fn)
    print(tyr_ata_sto)