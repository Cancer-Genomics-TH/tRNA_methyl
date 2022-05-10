#!/scratch/dkedra/soft/progs/pypy3.9_current/bin/pypy3

import sys
import gzip

input_fn = sys.argv[1]

def parse_list(fqgrep_columns):
    read_name = fqgrep_columns[0]
    read_seq  = fqgrep_columns[-2]
    read_qual = fqgrep_columns[-1]
    match = fqgrep_columns[-3]
    
    if match == "*":
        out_seq  = read_seq
        out_qual = read_qual
        #print(f"no_match: {fqgrep_columns}")
    else:
        start = int(sl[-5])
        end   = int(sl[-4])
        #print(start, end)
        seq_mask  = "N" * (end-start)
        qual_mask = "B" * (end-start)
        out_seq =  f"{read_seq[:start]}{seq_mask}{read_seq[end:]}"
        out_qual = f"{read_qual[:start]}{qual_mask}{read_qual[end:]}"
        #print(out_seq) 
        #print(f"got_match: {fqgrep_columns}")
    print(f"@{read_name}")
    print(out_seq)
    print("+")
    print(out_qual)


with gzip.open(input_fn, mode="rt") as fh:
    counter = 0
    for line in fh:
        if counter == 0:
            pass
        else:
            line = line.strip()
            sl = line.split()
            parse_list(sl)
            #print(sl)
        counter += 1
