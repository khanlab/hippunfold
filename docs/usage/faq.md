# Frequently asked questions

1. [](run-inference-mem)
2. [](no-input-images)
3. [](container-size)
4. [](model-files)


(run-inference-mem)=
## Why is the workflow stopping at the run_inference step?

If you are getting an error in the run_inference step, e.g. as follows:
```
[Thu Nov 10 02:11:20 2022]
Finished job 65.
18 of 193 steps (9%) done
Select jobs to execute...

[Thu Nov 10 02:11:20 2022]
rule run_inference:
    input: work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-preproc_T1w.nii.gz, /opt/hippunfold_cache/trained_model.3d_fullres.Task101_hcp1200_T1w.nnUNetTrainerV2.model_best.tar
    output: work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-nnunet_dseg.nii.gz
    log: logs/sub-1425/sub-1425_hemi-R_space-corobl_nnunet.txt
    jobid: 64
    reason: Missing output files: work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-nnunet_dseg.nii.gz; Input files updated by another job: work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-preproc_T1w.nii.gz
    wildcards: subject=1425, hemi=R
    resources: tmpdir=/tmp, gpus=0, mem_mb=16000, time=60

mkdir -p tempmodel tempimg templbl && cp work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-preproc_T1w.nii.gz tempimg/temp_0000.nii.gz && tar -xf /opt/hippunfold_cache/trained_model.3d_fullres.Task101_hcp1200_T1w.nnUNetTrainerV2.model_best.tar -C tempmodel && export RESULTS_FOLDER=tempmodel && export nnUNet_n_proc_DA=1 && nnUNet_predict -i tempimg -o templbl -t Task101_hcp1200_T1w -chk model_best --disable_tta &> logs/sub-1425/sub-1425_hemi-R_space-corobl_nnunet.txt && cp templbl/temp.nii.gz work/sub-1425/anat/sub-1425_hemi-R_space-corobl_desc-nnunet_dseg.nii.gz
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2022-11-10T020645.651622.snakemake.log
```
it is likely that you do not have enough memory available on your system. You need to have at least 8GB of memory on your system. If you are running Docker on Windows/Mac or another 
virtual machine (e.g. VirtualBox) you will need to increase the amount of memory dedicated to the virtual machine.

(no-input-images)=
## Why do I get the error, `No input images found for T1w`, or `No input images found for T2w`

The workflow is unable to find any input files to run HippUnfold.

This can happen if:
 - Singularity or docker cannot access your input directory. For Singularity, ensure your [Singularity options](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) are appropriate, in particular `SINGULARITY_BINDPATH`. For docker, ensure you are mounting the correct directory with the `-v` flag described in the [Getting started](https://hippunfold.readthedocs.io/en/latest/getting_started/docker.html) section. 
 - HippUnfold does not recognize your BIDS-formatted input images. This can occur if, for example, T1w images are labelled with the suffix `_t1w.nii.gz` instead of `_T1w.nii.gz` as per [BIDS specifications](https://bids.neuroimaging.io/specification.html). HippUnfold makes use of [PyBIDS](https://github.com/bids-standard/pybids) to parse the dataset, so we suggest you use the [BIDS Validator](https://bids-standard.github.io/bids-validator/) to ensure your dataset has no errors. Note: You can override BIDS parsing and use custom filenames with the `--path-*` option as described in the [](../usage/useful_options.md#parsing-non-bids-datasets-with-custom-paths) section. 

(container-size)=
## Why is the HippUnfold Docker/Singularity/Apptainer container so large?

In addition to some large software dependencies, the container has historically included U-net models for all the possible modalities we trained, each model taking up 2-4GB. 
We have addressed this issue in versions >= 1.3.0, by updating the workflow to download models on the fly (when they have not been previously downloaded), and not including any 
models in the container itself. This drops the container size significantly (<4GB compressed).

(model-files)=
## Why do I end up with large files in `~/.cache/hippunfold` after running HippUnfold? 

This folder is where the nnU-net model parameters are stored by default. You can override the location with the `HIPPUNFOLD_CACHE_DIR` environment variable. See [](../contributing/contributing.md#deep-learning-nnu-net-model-files) for more details.



