import nibabel as nib

ref_nib = nib.load(snakemake.input.ref_metric)
in_nib = nib.load(snakemake.params.in_metric)

vol = in_nib.get_fdata()
# vol = vol.reshape((vol.shape[0],vol.shape[1],1))

print(in_nib.affine)
print(ref_nib.affine)
print(vol.shape)

nib.Nifti1Image(vol, header=ref_nib.header, affine=ref_nib.affine).to_filename(
    snakemake.output.metric
)
