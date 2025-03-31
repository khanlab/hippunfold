import pyvista as pv
import numpy as np
import nibabel as nib
import nibabel.gifti as gifti
from nibabel.gifti.gifti import intent_codes
from collections import defaultdict
from lib.utils import setup_logger
from lib.surface import read_surface_from_gifti, find_boundary_vertices, write_label_gii
import sys

# Setup logger
log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)
sys.stdout = open(log_file, "a")
sys.stderr = open(log_file, "a")

logger.info("Loading surface from GIFTI...")
surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
logger.info(f"Surface loaded: {surface.n_points} vertices, {surface.n_faces} faces.")


logger.info("Find boundary vertices")
boundary_indices = find_boundary_vertices(surface)

boundary_scalars = np.zeros(surface.n_points, dtype=np.int32)  # Default is 0
boundary_scalars[boundary_indices] = 1  # Set boundary vertices to 1
logger.info(
    f"Boundary scalar array created. {np.sum(boundary_scalars)} boundary vertices marked."
)


# Find the largest connected component within this sub-mesh
logger.info("Applying largest connected components")

# Extract points that are within the boundary scalars
sub_mesh = pv.PolyData(surface.points, surface.faces).extract_points(
    boundary_scalars.astype(bool), adjacent_cells=True
)

# Compute connectivity to find the largest connected component
connected_sub_mesh = sub_mesh.connectivity("largest")

# Get indices of the largest component in the sub-mesh
largest_component_mask = (
    connected_sub_mesh.point_data["RegionId"] == 0
)  # Largest component has RegionId 0
largest_component_indices = connected_sub_mesh.point_data["vtkOriginalPointIds"][
    largest_component_mask
]

# Create an array for all points in the original surface
boundary_scalars = np.zeros(surface.n_points, dtype=np.int32)

# Keep only the largest component
boundary_scalars[largest_component_indices] = 1


logger.info("Saving GIFTI label file...")

label_dict = {
    "Background": {"key": 0, "red": 1.0, "green": 1.0, "blue": 1.0, "alpha": 0.0},
    "Boundary": {"key": 1, "red": 1.0, "green": 0.0, "blue": 0.0, "alpha": 1.0},
}

write_label_gii(
    boundary_scalars,
    snakemake.output.label_gii,
    label_dict=label_dict,
    metadata=metadata,
)
logger.info(f"GIFTI label file saved as '{snakemake.output.label_gii}'.")
