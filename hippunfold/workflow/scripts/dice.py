import nibabel as nib
import numpy as np

r = nib.load(snakemake.params.ref)
ref_mask = r.get_fdata()
n = nib.load(snakemake.input.res_mask)
tissuepred_rois = n.get_fdata()

tissuepred_bin = np.zeros(tissuepred_rois.shape)
lbls = snakemake.params.hipp_lbls
for l in lbls:
    tissuepred_bin[tissuepred_rois == l] = 1

dice = np.sum(ref_mask[tissuepred_bin == 1]) * 2.0 / (np.sum(ref_mask) + np.sum(tissuepred_bin))

# write to txt file
with open(snakemake.output.dice, "w") as f:
    f.write(str(dice))
