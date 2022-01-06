import numpy as np
import nibabel as nib

logfile = open(snakemake.log[0], "w")
print(f"start", file=logfile, flush=True)

## first load labelmap and binarize one edge that includes the inner surf and one edge the outer
lbl_nib = nib.load(snakemake.input[0])
lbl = lbl_nib.get_fdata()
print(f"labelmap loaded", file=logfile, flush=True)

source = np.zeros(lbl.shape)
for i in snakemake.params.src_labels:
    source[lbl == i] = 1
bin_nib = nib.Nifti1Image(source, lbl_nib.affine, lbl_nib.header)
nib.save(bin_nib, snakemake.output.innerbin)
print(f"inner edge binarized", file=logfile, flush=True)

sink = np.zeros(lbl.shape)
sink[lbl > 0] = 1
bin_nib = nib.Nifti1Image(sink, lbl_nib.affine, lbl_nib.header)
nib.save(bin_nib, snakemake.output.outerbin)
print(f"outer edge binarized", file=logfile, flush=True)
