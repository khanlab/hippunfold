import nibabel as nib

surf = nib.load(snakemake.input.surf_gii)
vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
vertices[:, 2] = snakemake.params.z_level
nib.save(surf, snakemake.output.surf_gii)
