# hippunfold
BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and subfield segmentation)

## Requirements:

Docker (Mac/Windows/Linux) or Singularity (Linux)

BIDS dataset with T1w and T2w images. Highly recommend using 0.8mm isotropic or higher 3D T2w TSE. 
Note: T1w only is available too with `--modality T1w`, however, this is discouraged unless you have high resolution (~0.7mm or better) T1w data, and performance will likely be sub-optimal. This is currently being evaluated.
Note 2: dwi workflows are also available but currently experimental

## How to run

### Docker:

Pull the container:

    docker pull khanlab/hippunfold:latest
do a dry run, printing the command at each step:

    docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -np 
run it with maximum number of cores:

    docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -p --cores all

### Singularity:

Pull the container:

    singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest
do a dry run, printing the command at each step:

    singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 
run it with maximum number of cores:

    singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -p --cores all

### Snakemake & Snakebids

This BIDS app makes use of [snakebids](https://github.com/akhanf/snakebids) to create a BIDS App from a [Snakemake](https://snakemake.github.io/) workflow. Thus, it can take any command-line arguments that snakemake can take, please see [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html), e.g. the `--dry-run, -n` and `--cores,-j` options are snakemake options. You can use a snakemake execution profile to submit jobs on a cluster as well (e.g. one for Compute Canada is [here](https://github.com/khanlab/cc-slurm))

## Setting up a dev environment:
```

Here are some instructions to get your python environment set-up on graham to run hippunfold:

1. create a virtualenv and activate it:
```
mkdir $SCRATCH/hippdev
cd $SCRATCH/hippdev
module load python/3
virtualenv venv
source venv/bin/activate
```
2. clone the source repos (so you can make/pull changes easily, or change to branch):
```
git clone --recursive http://github.com/khanlab/hippunfold
```

3. install snakebids using pip, with the -e option (for development mode):
```
pip install -e ./hippunfold
```

Now hippunfold will be installed for you and can run with:

hippunfold  <args here> 

Any containers used are included in the hippunfold workflow, and if in khanlab group on graham, will already be good to go..  If you log out, you just need to re-activate the virtualenv to start again 

If you ever want the latest code, can just pull it:
```
cd hippunfold
git pull
```

or if you need a branch, can: `git checkout <name of branch>`
```

### running on graham:

In an interactive job (for testing):

    regularInteractive -n 8
    hippunfold bids_dir out_dir participant --participant_label CC110037 -j 8

Submitting a job (for larger cores, more subjects), still single job, but snakemake will parallelize over the 32 cores:

    regularSubmit -j Fat hippunfold bids_dir out_dir participant  -j 32
 
Scaling up to ~hundred subjects (needs cc-slurm snakemake profile installed), submits 1 16core job per subject:

    hippunfold bids_dir out_dir participant  --profile cc-slurm

Scaling up to even more subjects (uses group-components to bundle multiple subjects in each job), 1 32core job for N subjects (e.g. 10):

    hippunfold bids_dir out_dir participant  --profile cc-slurm --group-components subj=10


[![](https://images.microbadger.com/badges/version/khanlab/hippunfold.svg)](https://microbadger.com/images/khanlab/hippunfold "Get your own version badge on microbadger.com")
[![Documentation Status](https://readthedocs.org/projects/hippunfold/badge/?version=latest)](https://hippunfold.readthedocs.io/en/latest/?badge=latest)
[![CircleCI](https://circleci.com/gh/khanlab/hippunfold.svg?style=svg)](https://circleci.com/gh/khanlab/hippunfold)

Documentation can be found [here.](https://hippunfold.readthedocs.io/en/latest/)
