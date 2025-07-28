import pyvista as pv
import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree
from lib.utils import setup_logger
from lib.surface import (
    write_surface_to_gifti,
    apply_affine_transform,
    remove_nan_vertices,
)
from copy import deepcopy

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

# Load data
coords_img = nib.load(snakemake.input.coords)
coords = coords_img.get_fdata()
nan_mask = nib.load(snakemake.input.nan_mask).get_fdata()
sink_mask = nib.load(snakemake.input.sink_mask).get_fdata()
src_mask = nib.load(snakemake.input.src_mask).get_fdata()
affine = coords_img.affine

# Prepare grid and apply transformations
grid = pv.ImageData(
    dimensions=np.array(coords.shape) + 1, spacing=(1, 1, 1), origin=(0, 0, 0)
)
coords[nan_mask == 1] = np.nan
coords[src_mask == 1] = -0.1
coords[sink_mask == 1] = 1.1

grid.cell_data["values"] = coords.flatten(order="F")
grid = grid.cells_to_points("values")

# Generate isosurface
logger.info("Generating isosurface")
surface = grid.contour(
    [snakemake.params.threshold], method="contour", compute_scalars=True
)
logger.info(surface)

# remove the nan-valued vertices - this isn't handled in PolyData.clean()
logger.info("Removing nan-valued vertices")
surface = remove_nan_vertices(surface)
logger.info(surface)

logger.info("Cleaning surface")
surface = surface.clean()
logger.info(surface)

logger.info("removing vertices not near coords")
# tag vertices not in coords
V = np.floor(surface.points).astype(int)
coord_at_V = coords[V[:, 0], V[:, 1], V[:, 2]]
valid_mask = np.logical_and(coord_at_V > 0, coord_at_V < 1)
logger.info(f"found {np.sum(valid_mask)} valid vertices")
# largest connected component on invalid vertices
badVmesh = deepcopy(surface)
badVmesh.points[valid_mask] = np.nan
badVmesh = remove_nan_vertices(badVmesh)
badVmesh = badVmesh.extract_largest()
logger.info(f"NaNing {badVmesh.n_points} connected invalid vertices")
# Use KDTree to map badVmesh.points back to surface.points
tree = cKDTree(surface.points)
_, inds = tree.query(badVmesh.points, k=1)
surface.points[inds] = np.nan
# remove the nan-valued vertices
logger.info("Removing nan-valued vertices")
surface = remove_nan_vertices(surface)
logger.info(surface)

logger.info("Extracting largest connected component")
surface = surface.extract_largest()
logger.info(surface)

logger.info(f"Filling holes up to radius {snakemake.params.hole_fill_radius}")
surface = surface.fill_holes(snakemake.params.hole_fill_radius)
logger.info(surface)

# reduce # of vertices with decimation
logger.info(f"Decimating surface with {snakemake.params.decimate_opts}")
surface = surface.decimate(snakemake.params.decimate_opts)
logger.info(surface)


logger.info(f"final surface clean to remove overlapping vertices, etc.")
surface = surface.clean()

# apply affine to go from matrix to image space
surface = surface.transform(affine, inplace=False)

# Save the final mesh
write_surface_to_gifti(surface, snakemake.output.surf_gii)
