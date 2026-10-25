import nibabel as nib
import numpy as np
from lib.qc_srcsink import compute_qc
from lib.surface import (
    assign_boundary_terminals,
    read_surface_from_gifti,
    write_label_gii,
)
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


nmin = snakemake.params.min_terminal_vertices
ap_length_penalty = snakemake.params.ap_length_penalty

edges = nib.load(snakemake.input.edges).agg_data()
ap_src = nib.load(snakemake.input.ap_src).agg_data()
ap_sink = nib.load(snakemake.input.ap_sink).agg_data()
pd_src = nib.load(snakemake.input.pd_src).agg_data()
pd_sink = nib.load(snakemake.input.pd_sink).agg_data()

logger.info("Loading surface from GIFTI...")
# structure metadata for the outputs comes from the surface
surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
logger.info(f"Surface loaded: {surface.n_points} vertices, {surface.n_faces} faces.")

logger.info("Assigning labels apsrc, apsink, pdsrc, pdsink")

boundary_vertices = np.where(edges == 1)[0]
distances = np.vstack(
    (
        ap_src[boundary_vertices],
        ap_sink[boundary_vertices],
        pd_src[boundary_vertices],
        pd_sink[boundary_vertices],
    )
).T
labels, info = assign_boundary_terminals(
    distances,
    boundary_vertices,
    surface,
    nmin,
    ap_length_penalty=ap_length_penalty,
    logger=logger,
)

label_counts = {k: int(np.sum(labels == k)) for k in range(4)}
logger.info(["Boundary vertices per label:", label_counts])

# QC metrics and checks, verified (and visualized) by downstream rules
qc = compute_qc(
    distances,
    boundary_vertices,
    labels,
    info,
    hemi=snakemake.wildcards.hemi,
    entities={"hemi": snakemake.wildcards.hemi, "label": snakemake.wildcards.label},
    thresholds=snakemake.params.qc_thresholds,
)
for check in qc.failed():
    logger.warning(f"QC {check['severity']}: {check['name']}: {check['message']}")
qc.save(snakemake.output.qc)


ap_srcsink = np.zeros((len(edges)), dtype=np.int32)
idx_edges = np.where(edges == 1)[0]
ap_srcsink[idx_edges[labels == 0]] = 1
ap_srcsink[idx_edges[labels == 1]] = 2

pd_srcsink = np.zeros((len(edges)), dtype=np.int32)
pd_srcsink[idx_edges[labels == 2]] = 1
pd_srcsink[idx_edges[labels == 3]] = 2


# save to gifti:
label_dict = {
    "Background": {"key": 0, "red": 1.0, "green": 1.0, "blue": 1.0, "alpha": 0.0},
    "Source": {"key": 1, "red": 1.0, "green": 0.0, "blue": 0.0, "alpha": 1.0},
    "Sink": {"key": 2, "red": 0.0, "green": 0.0, "blue": 1.0, "alpha": 1.0},
}

write_label_gii(
    ap_srcsink, snakemake.output.ap, label_dict=label_dict, metadata=metadata
)
write_label_gii(
    pd_srcsink, snakemake.output.pd, label_dict=label_dict, metadata=metadata
)
