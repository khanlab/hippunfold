# Specialized scans

This tutorial will cover how HippUnfold can be applied to non-standard
data including ex-vivo scans, super-high resolution data (eg. \<0.3mm
isotropic), non-MRI 3D imaging data, or scans where a corresponding
whole-brain T1w image is not available.

We will show how the available flags can be adapted for these use-cases
with several worked examples.

## Case 1: super high resolution

In this example, we have only a limited field of view covering the
hippocampus, and the resolution and contrast do not closely match the
training data of HippUnfold (0.3-1.0mm isotropic T1w, T2w, or DWI data).
This could be ex-vivo MRI data, or it could even be 3D microscopy data
as in our recent [3D BigBrain
publication](https://www.sciencedirect.com/science/article/pii/S105381191930919X).
Thus we don\'t expect HippUnfold\'s inbuilt UNet to be successful in
segmenting hippocampal tissue before unfolding, and we do not want to
downsample our data to accommodate HippUnfold\'s usual UNet and
unfolding workflow in `space-corobl` (which consists of 0.3mm isotropic
resampling cropped coronally-oblique to the hippocampus).

This will require manual segmentation of hippocampal grey matter, SRLM,
and neighbouring structures, though in the future we hope to include
models trained with higher resolution data (and contrasts more common in
ex-vivo scanning). This should be done according to the protocol
outlined
[here](https://ars.els-cdn.com/content/image/1-s2.0-S1053811917309977-mmc1.pdf)
or, more recently, the video example
[here](https://www.youtube.com/watch?v=mUQJ2GUcnLU&t=1s). This manual
segmentation file should have the `_dseg` suffix.

Here is an example of what the input directory might look like:

    exvivo/
    └── sub-001/
        ├── sub-001_hemi-R_desc-hippo_T2w.nii.gz
        └── sub-001_hemi-R_desc-hippo_dseg.nii.gz

This can be unfolded with the command:

    hippunfold . PATH_TO_OUTPUT_DIR participant --modality cropseg \
    --path_cropseg exvivo/sub-{subject}/sub-{subject}_hemi-{hemi}_desc-hippo_dseg.nii.gz \
    --hemi R --skip_inject_template_labels

Explanation: `--modality cropseg` informs HippUnfold that the input
manual segmentation should not be resampled and UNet does not need to be
run. Because of a limitation in bids parsing for the `hemi`
entity, we need to use the generic path input,
`--path_cropseg` in this case, making sure we use the
`{subject}` and `{hemi}` wildcards in the
filename. Output files will be named with `space-corobl` because
HippUnfold is coded to effectively treat all files as already being in
this space. We need the `--hemi R` to prevent HippUnfold looking for
both hemispheres. Finally, because this segmentation was performed
manually on very high resolution data, we can optionally consider
skipping the template shape injection step with
`--skip_inject_template_labels`. Template shape injection can fix minor
errors in segmentation from UNet or from an imperfect manual rater, at
the cost of smoothing out some details of the hippocampus due to the
fact that it uses deformable registration with inherent smoothness
contraints.

Note that because we are not resampling to the CITI168 template or using
UNet, the T2w image in this example is effectively not being used at
all. Instead, the provided manual segmentation makes up the basis for
unfolding.

## Case 2: one ex-vivo hemisphere 
In this example, we have a single hemisphere that was scanned ex-vivo at a
nearly standard resolution and T2w contrast. Because the resolution and
contrast are similar to the HippUnfold training data, we expect UNet
will work and so we don\'t need to perform manual segemntation. However,
due to gross deformations and the missing hemisphere, we don\'t expect
this sample to register well to the standard CITI168 template. Thus, we
will need to manually resample the image to the CITI168 template prior
to running HippUnfold, focusing in particular on aligning the
hippocampus. Once done, we may have a directory like this:

    PATH_TO_EXVIVO_DIR/
    └── sub-001/
        ├── sub-001_hemi-R_desc-exvivo_T2w.nii.gz
        ├── sub-001_hemi-R_affine-exvivo-to-CITI168_xfm.txt
        └── sub-001_hemi-R_desc-exvivo_space-CITI168_T2w.nii.gz

Note that only the last file is needed for unfolding:

    hippunfold . PATH_TO_OUTPUT_DIR participant --output_spaces corobl --hemi R --no_reg_template \
    --path_T2w PATH_TO_EXVIVO_DIR/sub-001/sub-001_hemi-R_desc-exvivo_space-CITI168_T2w.nii.gz \
    --output_spaces corobl

Here we need to use `--path_T2w` to specify which input should be used,
and `--no_reg_template` to specify that it is already in
`space-CITI168`. In this case, we also specified
`--output_spaces corobl`. This is not needed, but is useful when we are
interested in only the hippocampus as `space-corobl` is higher
resolution and cropped more nicely around the hippocampus then the
original scan, making it a good space to perform subsequent analyses.
Alternatively, outputs can be transformed back to the original space
using the inverted transform
`sub-001_hemi-R_affine-exvivo-to-CITI168_xfm.txt`.

This same usage could also be applied in a standard MRI case where no
T1w image is available.
