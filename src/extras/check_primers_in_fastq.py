#!/usr/bin/env python3

"""
script to run ripgrep (or grep) on fastq files
to check the number of reads from several FASTQ files
containing exact matches to primers from the short RNA kit:

https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/NEXTFLEX/miRNA/5132-08_NEXTflex_Small_RNA-Seq_v3_Automation_Guide_v17-04.pdf
"""

import glob

import executor

# primers_fn  = sys.argv[1]
primers_fn = "primers_NEXTflex.tsv"
primers_dir = "/data/sdb2/salamanca/proj/sandra_trna_20220423/github/tRNA_methyl/src"

primers_list = []

with open(f"{primers_dir}/{primers_fn}") as f:
    for line in f:
        line = line.strip()
        primer_name, primer_seq = line.split()
        primers_list.append((primer_name, primer_seq))

print(primers_list)

table_header_1 = "| fastq_name | num_of_reads |"
table_header_2 = " | ".join(x[0] for x in primers_list)
# print(table_header_1)
# print(table_header_2)
table_header = f"{table_header_1} | {table_header_2} |"
print(table_header)
table_header_bottom = "|"
for column_num in range(0, len(primers_list) + 2):
    table_header_bottom = f"{table_header_bottom}---|"
print(table_header_bottom)


fastq_fn_pattern = (
    "/data/sdb2/salamanca/proj/sandra_trna_20220423/FASTQ_data/*.fastq.gz"
)

for fastq_fn in glob.glob(fastq_fn_pattern):
    fastq_root = fastq_fn.split("/")[-1]
    column_values = [fastq_root]

    run_pigz = f"pigz --stdout --decompress {fastq_fn} | wc -l "
    num_of_reads = int(executor.execute(run_pigz, capture=True)) / 4
    column_values.append(num_of_reads)

    for primer_tuple in primers_list:
        primer_name, primer_seq = primer_tuple
        # print(fastq_fn, primer_name, primer_seq)
        command_1 = f"rg -z -c {primer_seq} {fastq_fn}"
        try:
            rg_out = executor.execute(command_1, capture=True)
            num_of_matches = int(rg_out)
        except:
            num_of_matches = 0

        column_values.append(num_of_matches)
    print("XXX", column_values)
    body_str = " | ".join(f"{x}" for x in column_values)
    out_str = f"| {body_str} |"
    print(out_str)
