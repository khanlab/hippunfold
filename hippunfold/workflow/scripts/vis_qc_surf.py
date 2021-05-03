from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib
matplotlib.use('Agg')

# 3D surface
# requires `snakemake.input.surf`
fig = plotting.plot_surf(snakemake.input.surf,view='dorsal')
fig.savefig(snakemake.output.png)
