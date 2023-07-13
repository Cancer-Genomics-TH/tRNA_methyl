#!/usr/bin/env python3

# coding: utf-8

""" 
env:
/scratch/dkedra/.conda/envs/python_310/lib/python3.10/

working version

downloads tRNA aligments as HTML from gtRNAdb
converts them to Stockholm format

output: hs38_gtrnadb_algn_01.sto
""" 


import executor
import requests
import sys

from bs4 import BeautifulSoup




url = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-align.html"

out_stockholm_fn = "trnas_hs38_gtrnadb_fromhtml.sto"

align_header = "# STOCKHOLM 1.0"




def parse_seq_alignments(alg_list):
    counter = 0
    old_id_prefix = ""
    old_anticodon = ""
    for alig in alg_list:
        lines = alig.text.split("\n")
        for line in lines:
            if line.find(">>") == -1 and line.find("filtered") == -1:
                sl = line.split()
                if len(sl) != 0:
                    #print(line)
                    seq, seq_id = sl[:2]
                    if seq_id.startswith("tRNA"):
                        split_id = seq_id.split("-")
                        id_prefix = "-".join(split_id[:3])
                        try:
                            anticodon = split_id[2]
                        except:
                            break
                            # there is an error in the HTML before the tRNA_Sup

                        if id_prefix != old_id_prefix:
                            if counter != 0:
                                print("//")
                            counter += 1
                            print(align_header)
                            align_group_id = f"#=GF ID {id_prefix}"
                            print(align_group_id)

                        out_str = f"{seq_id}\t{seq}"
                        print(out_str)
                        old_id_prefix = id_prefix
                    else:
                        split_id = seq_id.split("_")
                        anticodon = split_id[1]
                        short_id = f"{split_id[0]}_{split_id[1]}_HUMAN"
                        out_str = f"{short_id}\t{seq}"
                        print(out_str)

    print("//")    

def output_stockholm(alig_list):
    
    saveout = sys.stdout
    output_fh = open(out_stockholm_fn, "w")
    sys.stdout = output_fh

    parse_seq_alignments(alig_list)

    sys.stdout = saveout
    output_fh.close()
    #remove repeated lines
    command_1 = f"""cat {out_stockholm_fn} | uniq > {out_stockholm_fn}.tmp"""
    executor.execute(command_1)
    command_2 = f"""mv {out_stockholm_fn}.tmp {out_stockholm_fn}"""
    executor.execute(command_2)


if __name__ == "__main__":
    r = requests.get(url)
    soup = BeautifulSoup(r.text, 'html.parser')

    # get all alignments from html
    alig_list = soup.find_all('div', attrs={'class':'seq_alignment'})
    
    output_stockholm(alig_list)
   
