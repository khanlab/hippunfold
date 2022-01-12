# BIDS whole-brain data

This tutorial will cover applications of HippUnfold to an entire
[BIDS-compliant dataset](https://bids.neuroimaging.io/), meaning that
the same scan types are expected for all subjects which will be
processed in parallel. A typical call might look like this:

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant 

Depending on the method you used for installation, you may require
additional arguments such as `--cores all` or `--use-singularity`, or
prefixing the command with `singularity run`. This will expect
`PATH_TO_BIDS_DIR` to contain something like the following:

    PATH_TO_BIDS_DIR/
    └── sub-001/
        └── anat/
            ├── sub-001_T1w.nii.gz
            └── sub-001_T2w.nii.gz
    └── sub-002/
    ...

The `--modality` flag is required to specify which input image type should be used and in most cases, T1w should be most robust (though other types are supported!).

The T1w image is used to register to a standardized
template (CITI168), making it possible to reorient, upsample, and crop
around the left and right hippocampi (this is referred to within
HippUnfold as `space-corobl`). Note that the T1w image should 
have a whole-brain field of view. 

More examples of possible BIDS-compliant datasets can be found in
[hippunfold/test\_data/](https://github.com/khanlab/hippunfold/tree/master/test_data).

## Different input modalities 

By default, HippUnfold expects the `PATH_TO_BIDS_DIR` to contain at least
one T1w file for segmenting intrahippocampal
structures like the SRLM. However, we have
also provided models trained with T1w, T2w, or DWI data, or, users can input
their own custom manual segmentations for unfolding, which can be
specified with the `--modality` flag. For example:

    hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --modality T2w

would work for a dataset with only T2w images, like this one:

    PATH_TO_BIDS_DIR/
    └── sub-001/
        └── anat/
            └── sub-001_T2w.nii.gz
    ...

Note that in this case, registration to a T2w CITI168 template will be performed with the input T2w image. In some cases it may be preferrable to use a T1w image for registration to the standard CITI168 template. A T1w image can be registered to both the input T2w and T1w CITI168 template with the `--t1_reg_template` flag. This is typically most robust as long as a full brain FOV T1w image is available. If this registering is still failing then it may be improved with the `--rigid-reg-template` flag.

Specifying a manual segmentation (eg. `--modality segT1w`)
expects to additionally find an input file with the suffix `_dseg` which should
contain labels following the protocol outlined
[here](https://ars.els-cdn.com/content/image/1-s2.0-S1053811917309977-mmc1.pdf).
More details are provided on using manual segmentations on the following
page.

## Non-BIDS datasets

Wildcards can be used to enumarate input files if the data are not in
BIDS format. For example:

    PATH_TO_nonBIDS_DIR/
    └── sub-001_T1w.nii.gz
    └── sub-001_T2SPACE.nii.gz
    └── sub-001_TSE.nii.gz
    └── sub-002_T1w.nii.gz
    ...

This directory doesn\'t separate subjects into different folders or
contain an `anat/` folder for structural images. However, we can still
specify what subjects and images to use with `wildcards`. T2SPACE and
TSE are both acquisitions that are sensitive to T2-weights, but
HippUnfold will not recognize them without the suffix `_T2w`. We can
thus use the `--path_T2w` flag to specify exactly which of these file(s)
to use as inputs:

    hippunfold - PATH_TO_OUTPUT_DIR participant \
    --path_T1w PATH_TO_nonBIDS_DIR/sub-001_T1w.nii.gz \
    --path_T2w PATH_TO_nonBIDS_DIR/sub-{subject}_T2SPACE.nii.gz

This will search for any any files following the naming scheme and fill
in `{subject}` IDs for any files it can. Alternatively, `{subject}` IDs
can be provided in a list with the `--participant_label` flag.


