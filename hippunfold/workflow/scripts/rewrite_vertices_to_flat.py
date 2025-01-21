import nibabel as nib

gii = nib.load(snakemake.input.surf_gii)
vertices = gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data

ap = nib.load(snakemake.input.coords_AP).darrays[0].data
pd = nib.load(snakemake.input.coords_PD).darrays[0].data

vertices[:, 0] = ap * -float(snakemake.params.vertspace["extent"][0]
)  - float(snakemake.params.vertspace["origin"][0])
vertices[:, 1] = pd * float(snakemake.params.vertspace["extent"][1]
)  - float(snakemake.params.vertspace["origin"][1])
vertices[:, 2] = snakemake.params.z_level

nib.save(gii, snakemake.output.surf_gii)
