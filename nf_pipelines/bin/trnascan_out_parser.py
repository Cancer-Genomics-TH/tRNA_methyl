#!/usr/bin/env python3

"""
parsing trnascan-SE out file to get tRNA notes/annotations

"""

import io

import executor
import polars as pl


def parse_out(out_fn):
    """
    takes the original tRNAScan-SE output
    returns dictionary trna_id:tags
    """
    col_names = [
        "chrom",
        "trna_num",
        "trna_start",
        "trna_end",
        "trna_type",
        "anticodon",
        "intr_start",
        "intr_end",
        "inf_score",
        "iso_CM",
        "iso_score",
        "note",
    ]

    tags = ["pseudo", "IDP", "trunc_start"]

    command_1 = f"""tail -n +4 {out_fn} | tr -d ' ' """

    df = pl.read_csv(
        io.StringIO(executor.execute(command_1, capture=True)),
        has_header=False,
        infer_schema_length=2000,
        separator="\t",
        new_columns=col_names,
    ).filter(pl.col("chrom") != "MT")

    db_df = (
        df.with_columns(
            [
                pl.map(["chrom", "trna_num"], lambda s: s[0] + ".trna" + s[1]).alias(
                    "trna_id"
                ),
                pl.when(pl.col("note").is_not_null())
                .then(
                    pl.col("note").str.split(by=":").arr.first().str.replace(",", "_")
                )
                .otherwise("")
                .alias("tags"),
                pl.when(pl.col("intr_start") > 0)
                .then("intron")
                .otherwise("")
                .alias("introns"),
            ]
        )
        .filter((pl.col("tags") != "") | (pl.col("introns") != ""))
        .with_columns(
            pl.concat_str(pl.col(["tags", "introns"]), "_")
            .str.replace(r"^_|_$", "")
            .alias("final_tags")  # .fill_null("")
        )
    )
    # debug
    # print(db_df.filter(pl.col("tags") != "").head(200))
    my_df = db_df.select(["trna_id", "final_tags"])
    trna_ids = list(my_df["trna_id"])
    tags = list(my_df["final_tags"])
    tags_dict = dict(zip(trna_ids, tags))
    return tags_dict
    #  debug
    # print(tags_dict)
    # my_df.write_csv("test_trna_db_02.tsv", separator="\t")
    # filter((pl.col("tags") !="") | (pl.col("introns") !=""))


if __name__ == "__main__":
    trna_out_fn = "./results/chr_trnas_out"
    tags_dict_debug = parse_out(trna_out_fn)
    print(tags_dict_debug)
