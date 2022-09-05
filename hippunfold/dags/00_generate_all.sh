#!/bin/bash

hippunfold ../../test_data/bids_T1w /tmp/out participant -np --modality T1w --rulegraph --keep-work | ./generate_dag_images.sh T1w_rulegraph
hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --t1-reg-template --rulegraph --keep-work | ./generate_dag_images.sh T2w_t1-reg-template_rulegraph
hippunfold ../../test_data/bids_singleT2w /tmp/out participant -np --modality T2w --rulegraph --keep-work | ./generate_dag_images.sh T2w_rulegraph
hippunfold ../../test_data/bids_hippb500 /tmp/out participant -np --modality hippb500 --rulegraph --keep-work | ./generate_dag_images.sh hippb500_rulegraph
