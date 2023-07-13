#!/usr/bin/env python3

"""
converts mappings to selected tRNA contigs from a bam file to a fasta file with counts
"""

import io
import sys

import executor
import polars as pl

# import xxhash

trna_bed = "trnas_mature_all.bed"


bam_fn = sys.argv[1]
fasta_fn = f"""{bam_fn.split('.')[0]}.trna_mature.reads.hexhash.fa"""

samtools_com = f"""samtools view \
--min-MQ 1 \
--regions-file trnas_mature_all.bed  {bam_fn} \
| hck -f=3,4,10
"""

# foo = executor.execute( samtools_com, capture=True)
# print(foo)


df = pl.read_csv(
    io.StringIO(executor.execute(samtools_com, capture=True)),
    has_header=False,
    sep="\t",
    new_columns=["trna_name", "start", "seq"],
)

df = df.with_columns(pl.concat_str(pl.all()).hash(10, 20, 30, 40).alias("hashed"))

print(df.head())

df_counts = (
    df.groupby(["trna_name", "start", "seq", "hashed"])
    .count()
    .filter(pl.col("count") > 9)
)
print(df_counts.head())


df_counts = df_counts.sort(["trna_name", "start", "seq", "hashed"])

print(df_counts.head())


fa_df = df_counts.with_columns(
    (
        ">"
        + pl.col("trna_name")
        + "__strt"
        + pl.col("start").cast(pl.Utf8)
        + "__count"
        + pl.col("count").cast(pl.Utf8)
        + "__"
        + pl.col("hashed").cast(pl.Utf8)
        + "\n"
        + pl.col("seq")
    ).alias("fasta")
).drop(["trna_name", "seq", "start", "count", "hashed"])
fa_df.write_csv(fasta_fn, has_header=False, sep="\t")

sed_quotes = f"""sed -i 's/"//g'  {fasta_fn}"""
executor.execute(sed_quotes, capture=False)
