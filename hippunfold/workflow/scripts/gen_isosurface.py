import pyvista as pv
import nibabel as nib
import numpy as np
from copy import deepcopy
import pygeodesic.geodesic as geodesic
from scipy.ndimage import binary_dilation
from lib.utils import setup_logger

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


def write_surface_to_gifti(mesh, out_surf_gii):
    """
    Writes a PyVista mesh to a GIFTI surface file.
    """
    points = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices

    points_darray = nib.gifti.GiftiDataArray(
        data=points, intent="NIFTI_INTENT_POINTSET", datatype="NIFTI_TYPE_FLOAT32"
    )

    tri_darray = nib.gifti.GiftiDataArray(
        data=faces, intent="NIFTI_INTENT_TRIANGLE", datatype="NIFTI_TYPE_INT32"
    )

    gifti = nib.GiftiImage()
    gifti.add_gifti_data_array(points_darray)
    gifti.add_gifti_data_array(tri_darray)
    gifti.to_filename(out_surf_gii)


def apply_affine_transform(mesh, affine_matrix, inverse=False):
    """
    Applies an affine transformation to a PyVista mesh.
    """
    if inverse:
        affine_matrix = np.linalg.inv(affine_matrix)

    transformed_points = (
        np.hstack([mesh.points, np.ones((mesh.n_points, 1))]) @ affine_matrix.T
    )[:, :3]
    transformed_mesh = pv.PolyData(transformed_points, mesh.faces)
    return transformed_mesh


def remove_nan_vertices(mesh):
    """
    Removes NaN vertices from a PyVista mesh and updates faces accordingly.
    """
    valid_mask = ~np.isnan(mesh.points).any(axis=1)
    new_indices = np.full(mesh.n_points, -1, dtype=int)
    new_indices[valid_mask] = np.arange(valid_mask.sum())

    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices
    valid_faces_mask = np.all(valid_mask[faces], axis=1)
    new_faces = new_indices[faces[valid_faces_mask]]
    new_faces_pv = np.hstack([np.full((new_faces.shape[0], 1), 3), new_faces])

    cleaned_mesh = pv.PolyData(mesh.points[valid_mask], new_faces_pv)
    return cleaned_mesh


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
if snakemake.params.clean_method == "cleanAK":

    logger.info("Using cleanAK method")
    src_sink_nan_mask = get_adjacent_voxels(sink_mask, src_mask)
    coords[src_sink_nan_mask == 1] = np.nan


grid.cell_data["values"] = coords.flatten(order="F")
grid = grid.cells_to_points("values")
tfm_grid = grid.transform(affine, inplace=False)

# Generate isosurface
logger.info("Generating isosurface")
surface = tfm_grid.contour(
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

"""
#still experimenting with this..

logger.info("Applying surface smoothing")  
surface = surface.smooth_taubin(#normalize_coordinates=True,
                                #boundary_smoothing=True,
                                #feature_smoothing=True,
                                n_iter=20,
                                pass_band=0.1)  #n_iter=30, pass_band=0.1)
logger.info(surface)
"""

logger.info(f"Filling holes up to radius {snakemake.params.hole_fill_radius}")
surface = surface.fill_holes(snakemake.params.hole_fill_radius)
logger.info(surface)

# reduce # of vertices with decimation
logger.info(f"Decimating surface with {snakemake.params.decimate_opts}")
surface = surface.decimate_pro(**snakemake.params.decimate_opts)
logger.info(surface)


# logger.info("applying surface smoothing") # -- still experimenting with smoothing before/after decimation
# surface = surface.smooth_taubin()
# logger.info(surface)


if snakemake.params.clean_method == "cleanJD":

    ## JD clean - instead of trimming surfaces with a nan mask, we
    # keep vertices that overlap with good coord values. We then apply
    # some surface-based morphological opening and closing to keep
    # vertices along holes in the dg
    logger.info("Using cleanJD method")

    # this is equivalent to wb_command -volume-to-surface-mapping -enclosing
    # apply inverse affine to surface to get back to matrix space
    xfm_surface = apply_affine_transform(surface, affine, inverse=True)

    V = np.round(xfm_surface.points).astype("int") - 1
    # sample coords
    coord_at_V = coords[V[:, 0], V[:, 1], V[:, 2]]

    # keep vertices that are in a nice coordinate range
    epsilon = snakemake.params.coords_epsilon
    good_v = np.where((coord_at_V < (1 - epsilon)) & (coord_at_V > epsilon))[0]

    geoalg = geodesic.PyGeodesicAlgorithmExact(
        surface.points, surface.faces.reshape((-1, 4))[:, 1:4]
    )
    maxdist, _ = geoalg.geodesicDistances(good_v, None)
    bad_v = np.where(maxdist > snakemake.params.morph_openclose_dist)[0]
    maxdist, _ = geoalg.geodesicDistances(bad_v, None)
    bad_v = np.where(maxdist < snakemake.params.morph_openclose_dist)[0]

    surface.points[bad_v, :] = np.nan

    logger.info("Removing nan-valued vertices")
    surface = remove_nan_vertices(surface)
    logger.info(surface)

    logger.info("Extracting largest connected component")
    surface = surface.extract_largest()
    logger.info(surface)


# Save the final mesh
write_surface_to_gifti(surface, snakemake.output.surf_gii)
