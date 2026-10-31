from lib.surface_viewer import SurfaceViewer
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

viewer = SurfaceViewer(snakemake.params.title, entities=snakemake.params.entities)
for r in snakemake.params.records:
    viewer.add(
        r["subject"],
        r["hemi"],
        r["label"],
        r["native"],
        r["unfold"],
        r["metrics"],
        subfields=r["subfields"],
    )
logger.info(f"{len(snakemake.params.records)} surfaces")
viewer.save(snakemake.output.html)
