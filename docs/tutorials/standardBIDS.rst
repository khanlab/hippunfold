BIDS whole-brain data
=====================

This tutorial will cover applications of HippUnfold to an entire `BIDS-compliant dataset <https://bids.neuroimaging.io/>`_, meaning that the same scan types are expected for all subjects which will be processed in parallel. A typical call might look like this::

  hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant 
  

Depending on the method you used for installation, you may require additional arguments such as ``--cores all`` or ``--use-singularity``, or prefixing the command with ``singularity run``. This will expect ``PATH_TO_BIDS_DIR`` to contain something like the following::

  PATH_TO_BIDS_DIR/
  └── sub-001/
      └── anat/
          ├── sub-001_T1w.nii.gz
          └── sub-001_T2w.nii.gz
  └── sub-002/
  ...
          
          
In this case, the T1w image is used only to register to a standardized template (CITI168), making it possible to reorient, upsample, and crop around the left and right hippocampi (this is referred to within HippUnfold as ``space-corobl``). Note that only the T1w image needs to have a whole-brain field of view. By default, both of these input images are coregistered and preprocessed, but this can be skipped with the flags ``--skip_coreg`` and ``--skip_preproc``, repsectively, if this was already run on the data in ``PATH_TO_BIDS_DIR``. 

More examples of possible BIDS-complaint datasets can be found in `hippunfold/test_data/ <https://github.com/khanlab/hippunfold/tree/master/test_data>`_.

Different input modalities
------------------
By default, HippUnfold expects the ``PATH_TO_BIDS_DIR`` to contain at least one T2w file on which performance for segmenting intrahippocampal structures like the SRLM has been shown to be optimal. However, we have also provided models trained with T1w or DWI data, or, users can input their own custom manual segmentations for unfolding, which can be specified with the ``--modality`` flag. For example::

  hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --modality T1w
  

would work for a dataset with only T1w images, like this one::

  PATH_TO_BIDS_DIR/
  └── sub-001/
      └── anat/
          └── sub-001_T1w.nii.gz
  ...

Note that specifying a manual segmentation (eg. ``--modality segT1w``) expects to additionally find a file with the suffix ``_dseg`` which should contain labels following the protocol outlined `here <https://ars.els-cdn.com/content/image/1-s2.0-S1053811917309977-mmc1.pdf>`_. More details are provided on using manual segmentations on the following page.

Non-BIDS datasets
------------------
Wildcards can be used to enumarate input files if the data are not in BIDS format. For example::

  PATH_TO_nonBIDS_DIR/
  └── sub-001_T1w.nii.gz
  └── sub-001_T2SPACE.nii.gz
  └── sub-001_TSE.nii.gz
  └── sub-002_T1w.nii.gz
  ...

This directory doesn't separate subjects into different folders or contain an ``anat/`` folder for structural images. However, we can still specify what subjects and images to use with ``wildcards``. T2SPACE and TSE are both acquisitions that are sensitive to T2-weights, but HippUnfold will not recognize them without the suffix ``_T2w``. We can thus use the ``--path_T2w`` flag to specify exactly which of these file(s) to use as inputs::

  hippunfold - PATH_TO_OUTPUT_DIR participant \
  --path_T1w PATH_TO_nonBIDS_DIR/sub-001_T1w.nii.gz \
  --path_T2w PATH_TO_nonBIDS_DIR/sub-{subject}_T2SPACE.nii.gz

This will search for any any files following the naming scheme and fill in ``{subject}`` IDs for any files it can. Alternatively, ``{subject}`` IDs can be provided in a list with the ``--participant_label`` flag.

No T1w images
------------------
It is difficult to automatically reorient and crop images appropriately without a whole-brain T1w image, which we recommend collecting as a part of any acquisition protocol when possible. However, if this is not possible, the T2w image may be used instead by specifying it as the T1w with the ``--path_T1w`` flag. This registration may fail, and is especially likely to fail if the T2w images are not whole-brain. This can sometimes be ameliorated with the ``--rigid_reg_template`` flag. 

Alternatively, if you do not have a standard T1w scan, consider manually orienting and/or cropping your data and following the examples outlined in the `Specialized scans tutorial <https://github.com/khanlab/hippunfold/blob/tutorials/docs/tutorials/specializedScans.rst>`_.
