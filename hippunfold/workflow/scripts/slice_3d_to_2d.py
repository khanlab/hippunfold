import nibabel as nib

img = nib.load(snakemake.input.img)
matrix = img.get_fdata()[:, :, 0]
nib.Nifti1Image(matrix, img.affine).to_filename(snakemake.output.img)