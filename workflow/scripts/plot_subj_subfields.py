import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import seaborn as sns

df = pd.read_table(snakemake.input.tsv)

subjdf = df.drop(columns=['subject','hemi','Cyst']).transpose()
subjdf.columns=['L','R']

sns_plot = sns.lineplot(data=subjdf)
sns_plot.set_title(str(snakemake.wildcards))
sns_plot.set_ylabel('Volume (mm^3)')
sns_plot.get_figure().savefig(snakemake.output.png)

