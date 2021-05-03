from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib
matplotlib.use('Agg')

# original plot (centered)
fig = plotting.plot_roi(roi_img=snakemake.input.seg, bg_img=snakemake.input.img, display_mode='ortho',view_type='continuous',alpha=0.5,dim=-1,draw_cross=False)
fig.savefig(snakemake.output.png1)

# move 2mm backward in each direction
fig2 = plotting.plot_roi(roi_img=snakemake.input.seg, bg_img=snakemake.input.img, cut_coords=[x-2 for x in fig.cut_coords], display_mode='ortho',view_type='continuous',alpha=0.5,dim=-1,draw_cross=False)
fig2.savefig(snakemake.output.png2)

# move 2mm forward in each direction
fig3 = plotting.plot_roi(roi_img=snakemake.input.seg, bg_img=snakemake.input.img, cut_coords=[x+2 for x in fig.cut_coords], display_mode='ortho',view_type='continuous',alpha=0.5,dim=-1,draw_cross=False)
fig3.savefig(snakemake.output.png3)

