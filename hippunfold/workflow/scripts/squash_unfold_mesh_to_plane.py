import nibabel as nib
import numpy as np


def squash_gii_to_z_plane(input_gii_path, output_gii_path, z_value):
    """
    Modify the z-coordinate of all vertices in a surface GIFTI file to a fixed value.

    Parameters:
    - input_gii_path (str): Path to the input GIFTI file.
    - output_gii_path (str): Path to save the modified GIFTI file.
    - z_value (float): The fixed z-coordinate to set for all vertices.

    Returns:
    - None
    """
    # Load the input GIFTI file
    gii = nib.load(input_gii_path)

    # Extract the coordinates from the data array
    coords_array = gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0]
    coords = coords_array.data

    # Modify the z-coordinate
    coords[:, 2] = z_value

    # Update the data in the GIFTI object
    coords_array.data = coords

    # Ensure all metadata, including CIFTI-specific metadata, is preserved
    metadata = gii.meta
    coords_array.meta = coords_array.meta  # Preserve metadata in the array
    coords_array.dims = coords.shape  # Update dimensions in metadata

    # Save the modified GIFTI file
    nib.save(gii, output_gii_path)


squash_gii_to_z_plane(
    snakemake.input.surf_gii,
    snakemake.output.surf_gii,
    z_value=snakemake.params.z_level,
)
