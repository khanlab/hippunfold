import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

import pandas as pd
import seaborn as sns

df = pd.read_table(snakemake.input.csv, sep=",")

# subjdf = df.drop(columns=["subject", "hemi", "Cyst"]).transpose()
# subjdf.columns = ["L", "R"]

sns_plot = sns.barplot(data=df)
sns_plot.set_title(str(snakemake.wildcards))
sns_plot.set_ylabel(snakemake.wildcards.metric)
sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
sns_plot.get_figure().savefig(snakemake.output.png)
