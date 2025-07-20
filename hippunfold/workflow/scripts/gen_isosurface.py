import pyvista as pv
import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation
from lib.utils import setup_logger
from pymeshfix._meshfix import PyTMesh
from lib.surface import (
    write_surface_to_gifti,
    apply_affine_transform,
    remove_nan_vertices,
)

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def get_adjacent_voxels(mask_a, mask_b):
    """
    Create a mask for voxels where label A is adjacent to label B.
    Parameters:
    - mask_a (np.ndarray): A 3D binary mask for label A.
    - mask_b (np.ndarray): A 3D binary mask for label B.

    Returns:
    - np.ndarray: A 3D mask where adjacent voxels for label A and label B are marked as True.
    """
    # Dilate each mask to identify neighboring regions
    dilated_a = binary_dilation(mask_a)
    dilated_b = binary_dilation(mask_b)

    # Find adjacency: voxels of A touching B and B touching A
    adjacency_mask = (dilated_a.astype("bool") & mask_b.astype("bool")) | (
        dilated_b.astype("bool") & mask_a.astype("bool")
    )

    return adjacency_mask


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

# we also need to use a nan mask for the voxels where src and sink meet directly
# (since this is another false boundary)..
src_sink_nan_mask = get_adjacent_voxels(sink_mask, src_mask)
coords[src_sink_nan_mask == 1] = np.nan

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
surface = surface.clean(point_merging=False)
logger.info(surface)

logger.info("Extracting largest connected component")
surface = surface.extract_largest()
logger.info(surface)

logger.info(f"Filling holes up to radius {snakemake.params.hole_fill_radius}")
surface = surface.fill_holes(snakemake.params.hole_fill_radius)
logger.info(surface)

logger.info(f"filling more holes with mfix")
MeshFix = PyTMesh()
MeshFix.load_array(surface.points, surface.faces.reshape(-1, 4)[:, 1:])
MeshFix.fill_small_boundaries(nbe=100, refine=True)
vert, faces = MeshFix.return_arrays()
triangles = np.empty((faces.shape[0], 4), dtype=faces.dtype)
triangles[:, -3:] = faces
triangles[:, 0] = 3
surface = pv.PolyData(vert, triangles)
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
