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

logger.info("Loading surface from GIFTI...")
surface = read_surface_from_gifti(snakemake.input.surf_gii)
logger.info(f"Surface loaded: {surface.n_points} vertices, {surface.n_faces} faces.")


logger.info("applying strong surface smoothing") 
#surface = surface.smooth_taubin(n_iter=400,boundary_smoothing=False,feature_smoothing=False,normalize_coordinates=True,
#                                pass_band=0.01,
#                                edge_angle=45,#default 15
#                                feature_angle=90, #default 45
#                                )
surface = surface.smooth(**snakemake.params.smoothing_params)
logger.info(surface)


logger.info("writing output surface gifti") 
write_surface_to_gifti(surface, snakemake.output.smoothed_surf_gii)

