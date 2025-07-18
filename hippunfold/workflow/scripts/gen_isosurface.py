import pyvista as pv
import nibabel as nib
import numpy as np
from lib.utils import setup_logger
from lib.surface import (
    write_surface_to_gifti,
    apply_affine_transform,
    remove_nan_vertices,
)

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def get_second_largest_component(surface):
    # Extract all connected components
    components = surface.connectivity(split_bodies=True)
    
    # Get the component ID array
    component_ids = components['RegionId']
    
    # Count the number of faces in each component
    region_sizes = {}
    for region_id in range(component_ids.max() + 1):
        region_mask = component_ids == region_id
        region_size = region_mask.sum()
        region_sizes[region_id] = region_size
    
    # Sort regions by size in descending order
    sorted_regions = sorted(region_sizes.items(), key=lambda x: x[1], reverse=True)
    
    if len(sorted_regions) < 2:
        raise ValueError("Surface does not contain at least two connected components.")
    
    # Extract the second largest region ID
    second_largest_region_id = sorted_regions[1][0]
    
    # Threshold to extract just this region
    second_largest_component = components.threshold(value=second_largest_region_id, scalars='RegionId', preference='cell')
    
    return second_largest_component


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

# keep second largest connected component
logger.info("Keeping second largest CC (first should be outer surface)")
surface = get_second_largest_component(surface)
logger.info(surface)

logger.info("Cleaning surface")
surface = surface.clean()
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
