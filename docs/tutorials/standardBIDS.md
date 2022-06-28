# BIDS whole-brain data

This tutorial will cover applications of HippUnfold to an entire
[BIDS-compliant dataset](https://bids.neuroimaging.io/), meaning that
the same scan types are expected for all subjects which will be
processed in parallel. A typical call might look like this:

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --modality T1w

Depending on the method you used for installation, you may require
additional arguments such as `--cores all` or `--use-singularity`, or
prefixing the command with `singularity run`. This will expect
`PATH_TO_BIDS_DIR` to contain something like the following:

    PATH_TO_BIDS_DIR/
    └── dataset_description.json
    └── sub-001/
        └── anat/
            ├── sub-001_T1w.nii.gz
            └── sub-001_T2w.nii.gz
    └── sub-002/
    ...

The `--modality` flag is **required** to specify which input image type should be used and in most cases, T1w should be most robust (though other types are supported!).

The T1w image is used to register to a standardized
template (CITI168), making it possible to reorient, upsample, and crop
around the left and right hippocampi (this is referred to within
HippUnfold as `space-corobl`). Note that the T1w image should 
have a whole-brain field of view. 

More examples of possible BIDS-compliant datasets can be found in
[hippunfold/test\_data/](https://github.com/khanlab/hippunfold/tree/master/test_data).

## Different input modalities 

If using `modality=T1w`, HippUnfold expects the `PATH_TO_BIDS_DIR` to contain at least
one T1w file for segmenting intrahippocampal
structures like the SRLM. However, we have
also provided models trained with other contrasts, such as T2w, or, users can input
their own custom manual segmentations for unfolding. These can be
specified with the `--modality` flag. For example:

    hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --modality T2w

would work for a dataset with only T2w images, like this one:

    PATH_TO_BIDS_DIR/
    └── dataset_description.json
    └── sub-001/
        └── anat/
            └── sub-001_T2w.nii.gz
    ...

Note that in this case, registration to a T2w CITI168 template will be performed with the input T2w image. In some cases it may be preferrable to use a T1w image for registration to the standard CITI168 template. A T1w image can be registered to both the input T2w and T1w CITI168 template with the `--t1-reg-template` flag. This is typically most robust as long as a full brain FOV T1w image is available. If this registering is still failing then it may be improved with the `--rigid-reg-template` flag.

Specifying a manual segmentation (eg. `--modality segT1w`)
expects to additionally find an input file with the suffix `_dseg` which should
contain labels following the protocol outlined
[here](https://ars.els-cdn.com/content/image/1-s2.0-S1053811917309977-mmc1.pdf).
More details are provided on using manual segmentations on the following
page.

## Selecting and excluding subjects to process

By default, hippunfold will run on **all** the subjects in a dataset. If you want to run only on a subset of subjects, you can use the `--participant_label` flag, e.g. adding:

    --participant-label 001 

would run only on `sub-001`. You can add additional subjects by listing additional arguments to this option, e.g.:

    --participant-label 001 002

runs for `sub-001` and `sub-001`.

Also, if you want to exclude a subject, you can use the `--exclude-participant-label` option.

## Parsing Non-BIDS datasets with custom paths

Custom paths can be used to parse input datasets if the data are not in
BIDS format, but still are uniquely identified by subject (or subject+session) identifiers. For example:

    PATH_TO_nonBIDS_DIR/
    └── s_001_T1w.nii.gz
    └── s_001_T2SPACE.nii.gz
    └── s_001_TSE.nii.gz
    └── s_002_T1w.nii.gz
    └── s_002_T2SPACE.nii.gz
    └── s_002_TSE.nii.gz
    ...

This directory doesn't separate subjects into different folders or
contain an `anat/` folder for structural images. However, we can still
specify what subjects and images to use with subject `wildcards`. This is done by using the `--path-{modality}` options to specify the absolute location of the `nii.gz` files. Note that here, T2SPACE and
TSE are both T2-weighted acquisitions, and can be captured by using the `--path-T2w` flag to specify exactly which of these file(s)
to use as inputs. For example, the following command:

    hippunfold - PATH_TO_OUTPUT_DIR participant \
    --modality T2w \
    --t1-reg-template \
    --path_T1w PATH_TO_nonBIDS_DIR/s_{subject}_T1w.nii.gz \
    --path_T2w PATH_TO_nonBIDS_DIR/s_{subject}_T2SPACE.nii.gz

will search for any files following the naming scheme and fill
in `{subject}` IDs for any files it finds, using the T1w and T2SPACE images for `T1w` and `T2w` inputs

### Prerequisities for using custom path parsing: 

Not all non-BIDS datasets can be parsed, and may still need some reformatting or renaming.

Specifically:
 - The subject (or subject/session) wildcard(s) can only contain letters or numbers, e.g. they cannot include underscores, hyphens, or spaces.
 - The subject (or subject/session) wildcard(s) must be the only unique identifiers in the filenames. 

For example, this datasets would be **ineligible**:

    PATH_TO_nonBIDS_DIR/
    └── s_2019-05-29_001_T1w.nii.gz
    └── s_2019-05-29_001_T2SPACE.nii.gz
    └── s_2019-05-29_001_TSE.nii.gz
    └── s_2018-02-24_002_T1w.nii.gz
    └── s_2018-02-24_002_T2SPACE.nii.gz
    └── s_2018-02-24_002_TSE.nii.gz
    ...
    
You would need to rename/symlink your images to remove the additional unique date identifiers, or integrate it into thes subject wildcard, ensuring only letters and numbers appear in the wildcard, e.g.:

    PATH_TO_nonBIDS_DIR/
    └── s_20190529s001_T1w.nii.gz
    └── s_20190529s001_T2SPACE.nii.gz
    └── s_20190529s001_TSE.nii.gz
    └── s_20180224s002_T1w.nii.gz
    └── s_20180224s002_T2SPACE.nii.gz
    └── s_20180224s002_TSE.nii.gz
    ...


