import pyvista as pv
import nibabel as nib
import numpy as np
from vtk import vtkNIFTIImageReader


def remove_nan_points_faces(vertices,faces):
    # Step 1: Identify valid vertices (no NaN values)
    nan_mask = np.isnan(vertices).any(axis=1)
    valid_vertices = ~nan_mask  # True for valid rows
    valid_indices = np.where(valid_vertices)[0]
    
    # Step 2: Create a mapping from old to new indices
    new_indices_map = -np.ones(vertices.shape[0], dtype=int)  # Default to -1 for invalid vertices
    new_indices_map[valid_indices] = np.arange(len(valid_indices))
    
    # Step 3: Update the faces array to remove references to invalid vertices
    # Replace old indices with new ones, and remove faces with invalid vertices
    new_faces = []
    for face in faces:
        # Map old indices to new ones
        mapped_face = new_indices_map[face]
        if np.all(mapped_face >= 0):  # Include only faces with all valid vertices
            new_faces.append(mapped_face)

    new_faces = np.array(new_faces)
    
    # Step 4: Remove invalid vertices from the array
    new_vertices = vertices[valid_vertices]

    return (new_vertices,new_faces)
    




# Load the NIfTI file
nifti_img = nib.load(snakemake.input.nii)
data = nifti_img.get_fdata()


# Load the mask
mask_img = nib.load(snakemake.input.mask)
mask = mask_img.get_fdata()


affine = nifti_img.affine

# Get voxel spacing from the header
voxel_spacing = nifti_img.header.get_zooms()[:3]

# Create a PyVista grid
grid = pv.ImageData()
grid.dimensions = np.array(data.shape) + 1  # Add 1 to include the boundary voxels
grid.spacing = (1, 1, 1)  # Use unit spacing and zero origin since we will apply affine
grid.origin = (0,0,0) 

# Add the scalar field
grid.cell_data["values"] = np.where(mask == 0, np.nan, data).flatten(order="F")
grid = grid.cells_to_points("values")

# apply the affine
tfm_grid = grid.transform(affine,inplace=False)  # Apply the rotation part of the affine


# the contour function produces the isosurface
surface = tfm_grid.contour([snakemake.params.threshold],method='contour')

surface = surface.decimate(float(snakemake.params.decimate_percent) / 100.0)

# faces from pyvista surface are formatted with number of verts each row
# reshape and remove the first col to get Nx3
faces = surface.faces
faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]

points = surface.points


#with nans in background we end up with nan vertices, we can remove
#these to end up with an open contour.. 
new_points,new_faces = remove_nan_points_faces(points,faces)


points_darray = nib.gifti.GiftiDataArray(
    data=new_points, intent="NIFTI_INTENT_POINTSET", datatype="NIFTI_TYPE_FLOAT32"
)

tri_darray = nib.gifti.GiftiDataArray(
    data=new_faces, intent="NIFTI_INTENT_TRIANGLE", datatype="NIFTI_TYPE_INT32"
)


gifti = nib.GiftiImage()
gifti.add_gifti_data_array(points_darray)
gifti.add_gifti_data_array(tri_darray)

gifti.to_filename(snakemake.output.surf_gii)
