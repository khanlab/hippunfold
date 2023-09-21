import pandas as pd

pd.concat([pd.read_table(in_tsv) for in_tsv in snakemake.input]).to_csv(
    snakemake.output.tsv, sep="\t", index=False
)
