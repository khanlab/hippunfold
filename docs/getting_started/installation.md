# Installation

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and
subfield segmentation)

## Requirements

-   Docker (Mac/Windows/Linux) or Singularity (Linux)
-   For those wishing to contribute or modify the code, `pip install` or `poetry install` are also available (Linux), but will still require singularity to handle some dependencies. See [Contributing to HippUnfold](https://hippunfold.readthedocs.io/en/latest/contributing/contributing.html).
-   GPU not required

### Notes:

-   Inputs to Hippunfold should typically be a BIDS dataset including T1w images. Higher-resolution data are preferred (\<= 0.8mm) but the pipeline will still work with 1mm T1w images. See [Tutorials](https://hippunfold.readthedocs.io/en/latest/tutorials/standardBIDS.html).
-   Other 3D imaging modalities (eg. ex-vivo MRI, 3D histology, etc.) can be used, but may require manual tissue segmentation as the current workflow relies on U-net segmentation trained only on common MRI modalities.

## Running with Docker

Pull the container:

    docker pull khanlab/hippunfold:latest

See HippUnfold usage docs:

    docker run -it --rm \
    khanlab/hippunfold:latest \
    -h

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
HippUnfold Docker container, and the fourth line contains the required
arguments for HippUnfold, after which you can additionally specify optional arguments. You may want to familiarize yourself with
[Docker options](https://docs.docker.com/engine/reference/run/), and an
overview of HippUnfold arguments is provided in the [Command line
interface](https://hippunfold.readthedocs.io/en/latest/usage/app_cli.html)
documentation section.

## Running with Singularity

Pull from dockerhub:

    singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest

See HippUnfold usage docs:

    singularity run -e khanlab_hippunfold_latest.sif -h

Do a dry run, printing the command at each step:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 

Run it with maximum number of cores:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -p --cores all

Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) to ensure this container can read and write to yout input and output directories, repsectively. 
