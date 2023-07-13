#!/usr/bin/env python3


"""
spld-Tyr-GTA-2-1

"""


import sys
import textwrap

from pyfaidx import Fasta

"""
    
            print(tmp_list, isotype)
            
            if spliced_flag == "intron":
                if isotype in spliced_trnas:
                    new_name = new_name + ".intron"
                fix_seq(record, new_name)
            elif spliced_flag == "spld":
                if isotype in spliced_trnas:
                    new_name = new_name + ".spld"
                    fix_seq(record, new_name)    
                else:
                    #new_name = new_name + ".discard"
                    pass
            else:
                print("wrong flag", spliced_flag)


                isotype_fasta = f"{isotype}.fa"
                print(isotype_fasta, isotype)

                saveout = sys.stdout
                with open(isotype_fasta, "w") as output_fh:
                    sys.stdout = output_fh

"""




def parse_gtrnadb_fa(fasta_in):
    spliced_trnas = "Arg-TCT Ile-TAT Leu-CAA Tyr-GTA Tyr-ATA".split()

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
        with open(isotype_fasta, "w") as output_fh:
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
     
     in_fa = sys.argv[1]
     parse_gtrnadb_fa(in_fa)



