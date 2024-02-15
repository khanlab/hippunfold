from nilearn import plotting
import matplotlib.pyplot as plt


import matplotlib

matplotlib.use("Agg")

display = plotting.plot_anat(snakemake.input.flo, display_mode="ortho", dim=-0.5)
display.add_contours(snakemake.params.ref, colors="r")
display.savefig(snakemake.output.png)
display.close()
