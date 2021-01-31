from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib
matplotlib.use('Agg')

fig = plotting.plot_roi(roi_img=snakemake.input.seg, bg_img=snakemake.input.img, display_mode='ortho',view_type='continuous',alpha=0.5,dim=-1,draw_cross=False)

fig.savefig(snakemake.output.png)

