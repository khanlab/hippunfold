#!/bin/bash


for graphtype in dag rulegraph
do

    hippunfold ../../test_data/bids_T1w /tmp/out participant -np --modality T1w --${graphtype} --keep-work | ./proc_subgraph.py out_${graphtype}/T1w
    hippunfold ../../test_data/bids_T1w /tmp/out participant -np --modality T1w --${graphtype} --keep-work --config autotop_labels=['hipp']  --hemi R -n | ./proc_subgraph.py out_${graphtype}/T1w_hemi-R_hipponly
    hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --${graphtype} --keep-work | ./proc_subgraph.py out_${graphtype}/T2w
    hippunfold ../../test_data/bids_multiT2w /tmp/out participant -np --modality T2w --${graphtype} --keep-work | ./proc_subgraph.py out_${graphtype}/T2w_multi
    hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --${graphtype} --t1-reg-template --keep-work | ./proc_subgraph.py out_${graphtype}/T2w_t1-reg-template
    hippunfold ../../test_data/bids_hippb500 /tmp/out participant -np --modality hippb500 --${graphtype} --keep-work | ./proc_subgraph.py out_${graphtype}/hippb500

done
