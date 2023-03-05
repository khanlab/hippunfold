from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib

matplotlib.use("Agg")

dim = (
    -0.5
)  # seems to be more reliable, dim=-1 was blacking out some images that had low dynamic range..


fig, (ax1, ax2, ax3) = plt.subplots(3, 1)


# original plot (centered)
display = plotting.plot_roi(
    axes=ax1,
    roi_img=snakemake.input.seg,
    bg_img=snakemake.input.img,
    display_mode="ortho",
    view_type="continuous",
    alpha=0.5,
    dim=dim,
    draw_cross=False,
)

# move 2mm backward in each direction
plotting.plot_roi(
    axes=ax2,
    roi_img=snakemake.input.seg,
    bg_img=snakemake.input.img,
    cut_coords=[x - 2 for x in display.cut_coords],
    display_mode="ortho",
    view_type="continuous",
    alpha=0.5,
    dim=dim,
    draw_cross=False,
)

# move 2mm forward in each direction
plotting.plot_roi(
    axes=ax3,
    roi_img=snakemake.input.seg,
    bg_img=snakemake.input.img,
    cut_coords=[x + 2 for x in display.cut_coords],
    display_mode="ortho",
    view_type="continuous",
    alpha=0.5,
    dim=dim,
    draw_cross=False,
)

# save as one figure
fig.savefig(snakemake.output.png)
