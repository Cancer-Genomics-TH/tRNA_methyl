#!/usr/bin/env pypy3


"""
env:
/scratch/dkedra/.conda/envs/python_310/

naming scheme: spld-Tyr-GTA-2-1

next step: files are used as input for mlocarna alignment

caveat: 
spld-Tyr-ATA is exception, because it has just one sequence 

"""


import sys
import textwrap

from pyfaidx import Fasta


spliced_trnas = "Arg-TCT Ile-TAT Leu-CAA Tyr-GTA Tyr-ATA".split()

def parse_gtrnadb_fa(fasta_in, dir_out):
    

    gtrnadb_fa = Fasta(fasta_in, as_raw=True, sequence_always_upper=True)

    isotypes_dict = {}

    old_isotype = ""
    for record in gtrnadb_fa:
        tmp_list = f"{record.name}".split("_tRNA-")[1].split("-")
        isotype = "-".join(tmp_list[:2])

        if isotype in spliced_trnas:
            if isotype != old_isotype:
                isotypes_dict[isotype] = []
            isotypes_dict[isotype].append(record.name)                
        else:
            pass
        old_isotype = isotype
    print(isotypes_dict)
    print(isotypes_dict.keys())
    for spliced_isotype in isotypes_dict.keys():
        isotype_fasta = f"{spliced_isotype}.fa"
        print(isotype_fasta, spliced_isotype)

        saveout = sys.stdout
        with open(f"./{dir_out}/tRNA_{isotype_fasta}", "w") as output_fh:
            sys.stdout = output_fh
            for seq_name in isotypes_dict[spliced_isotype]:
                record = gtrnadb_fa[seq_name]
                
                new_name = seq_name.replace("Homo_sapiens_tRNA", "spld")
                print(f">{new_name}")
                sequence = f"{record}".replace("U", "T")
                print(sequence)
        
            sys.stdout = saveout
            output_fh.close()
        

if __name__ == "__main__":
     
    in_fa = "hg38-mature-tRNAs.fa"
    dir_out = sys.argv[1]
    parse_gtrnadb_fa(in_fa, dir_out)



 
