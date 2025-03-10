#!/bin/bash
set -euxo pipefail
HIPPUNFOLD_CACHE_DIR=`pwd`/test_data/fake_models
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w 
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w --hemi R 
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w --hemi L 
./hippunfold/run.py test_data/bids_multiT2w test_out participant -np --modality T2w 
./hippunfold/run.py test_data/bids_T1w test_out participant -np --modality T1w 
./hippunfold/run.py test_data/bids_hippb500 test_out participant -np --modality hippb500 
./hippunfold/run.py test_data/bids_T1w_longitudinal test_out participant -np --modality T1w 
./hippunfold/run.py test_data/bids_singleT2w_longitudinal test_out participant -np --modality T2w 
./hippunfold/run.py test_data/bids_dsegtissue test_out participant -np --modality dsegtissue --path_dsegtissue test_data/bids_dsegtissue/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_desc-tissue_dseg.nii.gz
./hippunfold/run.py test_data/bids_dsegtissue_1hemi test_out participant -np --modality dsegtissue --hemi L
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w --t1_reg_template 
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w --output_space T1w 
./hippunfold/run.py test_data/bids_T1w test_out participant -np --modality T1w --use-template-seg 
./hippunfold/run.py test_data/bids_singleT2w test_out participant -np --modality T2w --generate-myelin-map
./hippunfold/run.py test_data/bids_dsegtissue test_out participant_create_template -np --modality dsegtissue --path_dsegtissue test_data/bids_dsegtissue/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_desc-tissue_dseg.nii.gz --path-dsegsubfields test_data/bids_dsegtissue/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_desc-subfields_dseg.nii.gz
