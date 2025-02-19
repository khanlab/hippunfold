import pyvista as pv
import nibabel as nib
import numpy as np
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


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


# Load data
mask_img = nib.load(snakemake.input.mask)
mask = mask_img.get_fdata()
affine = mask_img.affine

# Prepare grid and apply transformations
grid = pv.ImageData(
    dimensions=np.array(mask.shape) + 1, spacing=(1, 1, 1), origin=(0, 0, 0)
)


grid.cell_data["values"] = mask.flatten(order="F")
grid = grid.cells_to_points("values")
tfm_grid = grid.transform(affine, inplace=False)

# Generate isosurface
logger.info("Generating isosurface")
surface = tfm_grid.contour(
    [snakemake.params.threshold], method="contour", compute_scalars=True
)
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

# reduce # of vertices with decimation
logger.info(f"Decimating surface with {snakemake.params.decimate_opts}")
surface = surface.decimate_pro(**snakemake.params.decimate_opts)
logger.info(surface)



# Save the final mesh
write_surface_to_gifti(surface, snakemake.output.surf_gii)


