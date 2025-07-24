import pyvista as pv
import numpy as np
from lib.utils import setup_logger
from lib.surface import read_surface_from_gifti, find_boundary_vertices, write_label_gii

# Setup logger
log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

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

# Extract points that are within the boundary scalars
sub_mesh = pv.PolyData(surface.points, surface.faces).extract_points(
    boundary_scalars.astype(bool), adjacent_cells=True
)

# Compute connectivity to find all connected components in the sub-mesh
connected_sub_mesh = sub_mesh.connectivity()

# Extract RegionId (connected component labels)
region_ids = connected_sub_mesh.point_data["RegionId"]
num_components = region_ids.max() + 1  # RegionIds are 0-based

# Count the number of points in each component
component_sizes = np.bincount(region_ids)

logger.info(f"Found {num_components} connected components.")

# Identify the largest component
largest_component_id = component_sizes.argmax()
largest_component_mask = region_ids == largest_component_id

# Create final scalar array for label output
boundary_scalars = np.zeros(surface.n_points, dtype=np.int32)

# Compute hole radii for smaller components
hole_radii = []

for region_id, size in enumerate(component_sizes):
    logger.info(f"Component {region_id}: {size} vertices")

    if region_id == largest_component_id:
        continue  # Skip largest component

    # Mask of points in this region
    region_mask = region_ids == region_id
    point_ids = connected_sub_mesh.point_data["vtkOriginalPointIds"][region_mask]
    boundary_scalars[point_ids] = 2
    coords = surface.points[point_ids]

    # Estimate hole radius as max distance from centroid
    centroid = coords.mean(axis=0)
    dists = np.linalg.norm(coords - centroid, axis=1)
    radius = dists.max()
    hole_radii.append(radius)

    logger.info(f"  → Estimated radius of component {region_id}: {radius:.3f}")

# Map back to original surface point indices
largest_component_indices = connected_sub_mesh.point_data["vtkOriginalPointIds"][
    largest_component_mask
]

boundary_scalars[largest_component_indices] = 1

logger.info("Saving GIFTI label file...")

label_dict = {
    "Background": {"key": 0, "red": 1.0, "green": 1.0, "blue": 1.0, "alpha": 0.0},
    "Boundary": {"key": 1, "red": 0.0, "green": 1.0, "blue": 0.0, "alpha": 1.0},
    "Holes": {"key": 2, "red": 1.0, "green": 0.0, "blue": 0.0, "alpha": 1.0},
}

write_label_gii(
    boundary_scalars,
    snakemake.output.label_gii,
    label_dict=label_dict,
    metadata=metadata,
)
logger.info(f"GIFTI label file saved as '{snakemake.output.label_gii}'.")
