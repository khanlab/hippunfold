import nibabel as nib
from lib.surface_stats import summarize_surface_metrics, vertex_areas

points, faces = nib.load(snakemake.input.surf).agg_data(
    ("NIFTI_INTENT_POINTSET", "NIFTI_INTENT_TRIANGLE")
)
metrics = {
    metric: nib.load(metric_gii).agg_data()
    for metric, metric_gii in zip(snakemake.params.metrics, snakemake.input.metrics)
}

# subfield labels exist for the hipp only, so the dentate is summarized whole
if hasattr(snakemake.input, "subfields"):
    label_gii = nib.load(snakemake.input.subfields)
    parcels = label_gii.agg_data()
    names = label_gii.labeltable.get_labels_as_dict()
else:
    parcels, names = None, None

df = summarize_surface_metrics(
    metrics,
    vertex_areas(points, faces),
    parcels=parcels,
    names=names,
    name=snakemake.params.subfield,
)

# prefixed as in gen_volume_tsv, which also keeps them strings when read back
entities = {"subject": f"sub-{snakemake.wildcards.subject}"}
if snakemake.wildcards.get("session"):
    entities["session"] = f"ses-{snakemake.wildcards.session}"
entities["hemi"] = snakemake.wildcards.hemi
entities["label"] = snakemake.wildcards.label
for i, (entity, value) in enumerate(entities.items()):
    df.insert(i, entity, value)

df.to_csv(snakemake.output.tsv, sep="\t", index=False, na_rep="n/a")
