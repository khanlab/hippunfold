import nibabel as nib
import numpy as np
import pyvista as pv
from lib.qc_surface import compute_topology
from lib.surface import write_points_faces_to_gifti
from lib.surface_repair import fill_holes, repair_mesh
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def log_topology(stage, points, faces):
    t = compute_topology(points, faces)
    logger.info(
        f"{stage}: V={t['V']} F={t['F']} chi={t['chi']} components={t['components']} "
        f"boundary loops={t['boundary_loops']} genus={t['genus']:g} "
        f"misoriented={len(t['misoriented_faces'])} "
        f"non-manifold edges={len(t['nonmanifold_edges'])} "
        f"vertices={len(t['nonmanifold_vertices'])}"
    )


def repair(stage, points, faces):
    points, faces, counts = repair_mesh(points, faces)
    changes = {k: v for k, v in counts.items() if k != "orientable" and v}
    logger.info(f"repair ({stage}): {changes or 'no changes'}")
    if not counts["orientable"]:
        logger.warning(f"repair ({stage}): surface is not orientable")
    return points, faces


points, faces = nib.load(snakemake.input.surf_gii).agg_data()
points = np.asarray(points, dtype=float)
faces = np.asarray(faces, dtype=np.int64)
log_topology("isosurface", points, faces)

points, faces = repair("isosurface", points, faces)

# hole_fill_radius is in voxels of the coords image the isosurface came from
voxel_size = float(np.mean(nib.load(snakemake.input.coords).header.get_zooms()[:3]))
max_radius = snakemake.params.hole_fill_radius * voxel_size
logger.info(f"Filling holes up to radius {max_radius:.2f} mm")
points, faces, radii = fill_holes(points, faces, max_radius)
logger.info(
    f"filled {len(radii)} holes"
    + (f", radius (mm): {', '.join(f'{r:.2f}' for r in radii)}" if radii else "")
)
points, faces = repair("after filling holes", points, faces)
log_topology("after filling holes", points, faces)

logger.info(f"Decimating surface with {snakemake.params.decimate_opts}")
surface = pv.PolyData(points, np.c_[np.full(len(faces), 3), faces].ravel())
surface = surface.decimate(snakemake.params.decimate_opts).clean()
points = np.asarray(surface.points, dtype=float)
faces = surface.faces.reshape((-1, 4))[:, 1:].astype(np.int64)
log_topology("after decimation", points, faces)

points, faces = repair("after decimation", points, faces)
log_topology("repaired", points, faces)

write_points_faces_to_gifti(points, faces, snakemake.output.surf_gii)
