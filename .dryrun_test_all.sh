#!/bin/bash
set -euxo pipefail
HIPPUNFOLD_CACHE_DIR=`pwd`/test_data/fake_models
hippunfold test_data/bids_singleT2w test_out participant -np --modality T2w
hippunfold test_data/bids_singleT2w test_out participant -np --modality T2w --hemi R
hippunfold test_data/bids_singleT2w test_out participant -np --modality T2w --hemi L
hippunfold test_data/bids_multiT2w test_out participant -np --modality T2w
hippunfold test_data/bids_T1w test_out participant -np --modality T1w
hippunfold test_data/bids_hippb500 test_out participant -np --modality hippb500
hippunfold test_data/bids_T1w_longitudinal test_out participant -np --modality T1w
hippunfold test_data/bids_singleT2w_longitudinal test_out participant -np --modality T2w
hippunfold test_data/bids_segT2w test_out participant -np --modality segT2w
hippunfold . test_out participant -np --modality cropseg --path_cropseg test_data/data_cropseg/sub-{subject}_hemi-{hemi}_dseg.nii.gz
hippunfold . test_out participant -np --modality cropseg --path_cropseg test_data/data_cropseg_1hemi/sub-{subject}_hemi-{hemi}_dseg.nii.gz --hemi L
hippunfold test_data/bids_singleT2w test_out participant -np --modality T2w --t1_reg_template
hippunfold test_data/bids_singleT2w test_out participant -np --modality T2w --output_space T1w
hippunfold test_data/bids_T1w test_out participant -np --modality T1w --use-template-seg
