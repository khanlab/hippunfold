#!/bin/bash


for graphtype in dag rulegraph
do

#    hippunfold ../../test_data/bids_T1w /tmp/out participant -np --modality T1w --${graphtype}  | tail -n +3 | ./proc_subgraph.py out_${graphtype}/T1w
#    hippunfold ../../test_data/bids_T1w /tmp/out participant -np --modality T1w --${graphtype}  --config autotop_labels=['hipp']  --hemi R -n | tail -n +3 | ./proc_subgraph.py out_${graphtype}/T1w_hemi-R_hipponly
#    hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --${graphtype}  | tail -n +3 | ./proc_subgraph.py out_${graphtype}/T2w
#    hippunfold ../../test_data/bids_multiT2w /tmp/out participant -np --modality T2w --${graphtype}  | tail -n +3 | ./proc_subgraph.py out_${graphtype}/T2w_multi
#    hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --${graphtype} --t1-reg-template  | tail -n +3 | ./proc_subgraph.py out_${graphtype}/T2w_t1-reg-template
#    hippunfold ../../test_data/bids_hippb500 /tmp/out participant -np --modality hippb500 --${graphtype}  | tail -n +3 | ./proc_subgraph.py out_${graphtype}/hippb500
    hippunfold ../../test_data/bids_dsegtissue /tmp/out group_create_atlas -np --modality dsegtissue --path_dsegtissue ../../test_data/bids_dsegtissue/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_desc-tissue_dseg.nii.gz --path-dsegsubfields ../../test_data/bids_dsegtissue/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_desc-subfields_dseg.nii.gz --new-atlas-name mytestatlas --${graphtype} | tail -n +3 | ./proc_subgraph.py out_${graphtype}/create_atlas

done
