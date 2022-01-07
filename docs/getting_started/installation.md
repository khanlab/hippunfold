# Installation

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and
subfield segmentation)

## Requirements

-   Docker (Mac/Windows/Linux) or Singularity (Linux)
-   `pip install` is also available (Linux), but will still require
    singularity to handle some dependencies
-   For the T1w workflow: BIDS dataset with T1w images.
    Higher-resolution data are preferred (\<= 0.8mm) but the pipeline
    will still work with 1mm T1w images.
-   GPU not required

### Notes:

-   T1w and/or T2w input images are supported, and we recommend using
    sub-millimetric isotropic data for best performance.
-   Other 3D imaging modalities (eg. ex-vivo MRI, 3D histology, etc.)
    can be used, but may require manual tissue segmentation.

## Running with Docker

Pull the container:

    docker pull khanlab/hippunfold:latest

Do a dry run, printing the command at each step:

    docker run -it --rm \
    -v PATH_TO_BIDS_DIR:/bids:ro \
    -v PATH_TO_OUTPUT_DIR:/output \
    khanlab/hippunfold:latest \
    /bids /output participant -np 

Run it with maximum number of cores:

    docker run -it --rm \
    -v PATH_TO_BIDS_DIR:/bids:ro \
    -v PATH_TO_OUTPUT_DIR:/output \
    khanlab/hippunfold:latest \
    /bids /output participant -p --cores all

For those not familiar with Docker, the first three lines of this
example are generic Docker arguments to ensure it is run with the safest
options and has permission to access your input and output directories
(specified here in capital letters). The third line specifies the
HippUnfold Docker container, and the fourth line contains the requires
arguments for HippUnfold. You may want to familiarize yourself with
[Docker options](https://docs.docker.com/engine/reference/run/), and an
overview of HippUnfold arguments is provided in the [Command line
interface](https://hippunfold.readthedocs.io/en/latest/usage/app_cli.html)
documentation section.

## Running with Singularity

If you are using singularity, there are two different ways to run
`hippunfold`.

1.  Pull the fully-contained BIDS app container from
    `docker://khanlab/hippunfold:latest` and run it
2.  Clone the github repo, `pip install`, and run `hippunfold`, which
    will pull the required singularity containers as needed.

If you want to modify the workflow at all or run a development branch,
then option 2 is required.

Furthermore, if you want to run the workflow on a Cloud system that
Snakemake supports, or with a Cluster Execution profile, you should use
Option 2 (see snakemake documentation on [cluster
profiles](https://github.com/snakemake-profiles/doc) and [cloud
execution](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)).

### Option 1:

Pull the container:

    singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest

Do a dry run, printing the command at each step:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 

Run it with maximum number of cores:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -p --cores all

### Option 2:

This assumes you have created and activated a virtualenv or conda
environment first.

1.  Install hippunfold using pip:

        pip install hippunfold

2.  Run the following to download the U-net models:

        hippunfold_download_models

3.  Now hippunfold will be installed for you and can perform a dry-run
    with:

        hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np

4.  To run on all cores, and have snakemake pull any required
    containers, use:

        hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --cores all --use-singularity

## Additional instructions for Compute Canada

This section provides an example of how to set up a `pip installed` copy
of HippUnfold on CompateCanada\'s `graham` cluster.

### Setting up a dev environment on graham:

Here are some instructions to get your python environment set-up on
graham to run HippUnfold:

1.  Create a virtualenv and activate it:

        mkdir $SCRATCH/hippdev
        cd $SCRATCH/hippdev
        module load python/3
        virtualenv venv
        source venv/bin/activate

2.  Follow the steps above to `pip install` from github repository

### Running hippunfold jobs on graham:

Note that this requires
[neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) for
regularSubmit or regularInteractive wrappers, and the
[cc-slurm](https://github.com/khanlab/cc-slurm) snakemake profile for
graham cluster execution with slurm.

In an interactive job (for testing):

    regularInteractive -n 8
    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant \
    --participant_label 001 -j 8

Here, the last line is used to specify only one subject from a BIDS
directory presumeably containing many subjects.

Submitting a job (for larger cores, more subjects), still single job,
but snakemake will parallelize over the 32 cores:

    regularSubmit -j Fat \
    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -j 32

Scaling up to \~hundred subjects (needs cc-slurm snakemake profile
installed), submits 1 16core job per subject:

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant \
    --profile cc-slurm

Scaling up to even more subjects (uses group-components to bundle
multiple subjects in each job), 1 32core job for N subjects (e.g. 10):

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant \
    --profile cc-slurm --group-components subj=10
