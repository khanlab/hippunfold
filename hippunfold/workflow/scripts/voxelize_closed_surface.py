import pyvista as pv
import nibabel as nib
import numpy as np
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def read_surface_from_gifti(surf_gii):
    """Load a surface mesh from a GIFTI file."""
    surf = nib.load(surf_gii)
    vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
    faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")
    faces_pv = np.hstack([np.full((faces.shape[0], 1), 3), faces])  # PyVista format

    return pv.PolyData(vertices, faces_pv)


# Load data
ref_img = nib.load(snakemake.input.ref_nii)
ref_vol = ref_img.get_fdata()
affine = ref_img.affine

logger.info("Loading surface from GIFTI...")
surface = read_surface_from_gifti(snakemake.input.surf_gii)
logger.info(f"Surface loaded: {surface.n_points} vertices, {surface.n_faces} faces.")


# Prepare grid and apply transformations
grid = pv.ImageData(
    dimensions=np.array(ref_vol.shape) + 1, spacing=(1, 1, 1), origin=(0, 0, 0)
)


grid.cell_data["values"] = ref_vol.flatten(order="F")
grid = grid.cells_to_points("values")
tfm_grid = grid.transform(affine, inplace=False)

vox = tfm_grid.select_enclosed_points(surface,tolerance=1e-8)
nib.Nifti1Image(vox.point_data['SelectedPoints'].reshape(ref_vol.T.shape).T,affine=ref_img.affine, header=ref_img.header).to_filename(snakemake.output.mask_nii)

