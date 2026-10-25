import nibabel as nib
from lib.qc_surface import compute_qc
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

# an intermediate stage (e.g. the isosurface before repair) has only the
# midthickness, and no coords to check the unfolded mapping
surfaces = {
    name: nib.load(snakemake.input[name]).agg_data()
    for name in ("midthickness", "inner", "outer")
    if hasattr(snakemake.input, name)
}
if hasattr(snakemake.input, "coords_ap"):
    coords_ap = nib.load(snakemake.input.coords_ap).agg_data()
    coords_pd = nib.load(snakemake.input.coords_pd).agg_data()
else:
    coords_ap = coords_pd = None

entities = {k: snakemake.wildcards[k] for k in ("hemi", "label")}
result = compute_qc(
    surfaces, coords_ap, coords_pd, entities, snakemake.params.qc_thresholds
)
result.add_data("stage", snakemake.wildcards.get("stage", "final"))
for check in result.checks:
    logger.info(
        f"{'pass' if check['passed'] else check['severity']}: {check['name']}: "
        f"{check['message']}"
    )
result.save(snakemake.output.qc)
logger.info(f"{len(result.data['defects'])} defects, saved to {snakemake.output.qc}")
