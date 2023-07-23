import nibabel as nib
import numpy as np

gii = nib.load(snakemake.input.gii)
varr = gii.darrays[0].data
varr = np.tanh(varr)

nib.save(gii, snakemake.output.gii)
