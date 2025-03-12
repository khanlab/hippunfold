import nibabel as nib

ref_nib = nib.load(snakemake.input.ref_nii)
in_nib = nib.load(snakemake.input.nii)

vol = in_nib.get_fdata()


nib.Nifti1Image(vol, header=ref_nib.header, affine=ref_nib.affine).to_filename(
    snakemake.output.nii
)
