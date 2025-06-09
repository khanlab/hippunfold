import numpy as np
import nibabel as nib
from scipy.ndimage import binary_dilation, generate_binary_structure


def selective_dilation(
    input_nifti, output_nifti, src, sink, src_bg, sink_bg, structure_size=3
):
    # Load image
    img = nib.load(input_nifti)
    data = img.get_fdata().astype(np.int32)

    # Create structuring element
    structure = generate_binary_structure(3, 1)

    # Dilation: src -> src_bg
    src_mask = data == src
    src_bg_mask = data == src_bg
    src_dilated = binary_dilation(src_mask, structure=structure) & src_bg_mask
    data[src_dilated] = src

    # Dilation: sink -> sink_bg
    sink_mask = data == sink
    sink_bg_mask = data == sink_bg
    sink_dilated = binary_dilation(sink_mask, structure=structure) & sink_bg_mask
    data[sink_dilated] = sink

    # Save the modified image
    new_img = nib.Nifti1Image(data, img.affine, img.header)
    nib.save(new_img, output_nifti)

    print(f"Output saved to {output_nifti}")


input_nifti = snakemake.input.template_seg
output_nifti = snakemake.output.template_seg
src = snakemake.params.src_label
sink = snakemake.params.sink_label
src_bg = snakemake.params.src_bg
sink_bg = snakemake.params.sink_bg
structure_size = snakemake.params.struc_elem_size

selective_dilation(
    input_nifti, output_nifti, src, sink, src_bg, sink_bg, structure_size
)
