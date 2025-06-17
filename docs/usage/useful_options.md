# Running HippUnfold on your data

This section goes over the command-line options you will find
most useful when running HippUnfold on your dataset, along with 
describing some of the issues you might face.

Note: Please first refer to the simple example in the Installation 
section, which goes over running HippUnfold on a test dataset, and the
essential required options.

## Selecting the modality to use

The `--modality` option must be chosen when running HippUnfold, and it affects what 
U-net model will be used, and how the pre-processing will be performed on 
the images.

If you have sub-millimetric, isotropic, whole-brain T1w data, the `--modality T1w` option is recommended. 

If you T2w data, you can use the `--modality T2w` option, however, you may need to also
use the T1w data for template registration (`--t1-reg-template`), especially if you have a limited FOV.
This is typically most robust as long as a full brain FOV T1w image is available. If this registration
 is still failing then it may be improved with the `--rigid-reg-template` flag.

For protocols employing high-resolution, b-value 500, hippocampal diffusion-weighted imaging, 
the `--modality hippb500` option can be used, and does not require registration
to a template (providing your acquisition is axial and oblique to the hippocampus).

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


## Known limitations for BIDS parsing


HippUnfold uses snakebids, which makes use of pybids to parse 
a [BIDS-compliant dataset](https://bids.neuroimaging.io/). However,
because of the way Snakebids and Snakemake operate, one limitation is that 
the input files in your BIDS dataset need to be consistent in terms of 
what optional BIDS entities exist in them. We can use the acquisition (`acq`) 
entity as an example. HippUnfold should have no problem parsing the following dataset:

    PATH_TO_BIDS_DIR/
    └── dataset_description.json
    └── sub-001/
        └── anat/
            ├── sub-001_acq-mprage_T1w.nii.gz
    └── sub-002/
        └── anat/
            ├── sub-002_acq-spgr_T1w.nii.gz
    ...

as the path (with wildcards) will be interpreted as `sub-{subject}_acq-{acq}_T1w.nii.gz`.

However, the following dataset will raise an error:

    PATH_TO_BIDS_DIR/
    └── dataset_description.json
    └── sub-001/
        └── anat/
            ├── sub-001_acq-mprage_T1w.nii.gz
    └── sub-002/
        └── anat/
            ├── sub-002_T1w.nii.gz
    ...

because two distinct paths (with wildcards) would be found for T1w images:
```
sub-{subject}_acq-{acq}_T1w.nii.gz
``` 
and 
```
sub-{subject}_T1w.nii.gz
```

Similarly, you could not have some subjects with the `ses` identifier, 
and some subjects without it. 

There will soon be added functionality in snakebids to filter out extra files, 
but for now, if your dataset has these issues you will need to rename or remove extraneous files.

More examples of possible BIDS-compliant datasets can be found in
[hippunfold/test\_data/](https://github.com/khanlab/hippunfold/tree/master/test_data).


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

    hippunfold PATH_TO_nonBIDS_DIR PATH_TO_OUTPUT_DIR participant \
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


