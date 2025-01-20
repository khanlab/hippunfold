import nibabel as nib

gii = nib.load(snakemake.input.surf_gii)
vertices = gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data

ap = nib.load(snakemake.input.coords_AP).darrays[0].data
pd = nib.load(snakemake.input.coords_PD).darrays[0].data

vertices[:, 0] = ap
vertices[:, 1] = pd
vertices[:, 2] = 0  # TODO include IO?

nib.save(gii, snakemake.output.surf_gii)
