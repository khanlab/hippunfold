import matplotlib
from nilearn import plotting

matplotlib.use("Agg")

# 3D surface
# requires `snakemake.input.surf`
fig = plotting.plot_surf(snakemake.input.surf, view="dorsal")
fig.savefig(snakemake.output.png)
