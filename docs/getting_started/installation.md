# Installation

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and
subfield segmentation)

## Requirements

-   Docker (Intel Mac/Windows/Linux) or Singularity (Linux)
-   For those wishing to contribute or modify the code, `pip install` or `poetry install` are also available (Linux), but will still require singularity to handle some dependencies. See [Contributing to HippUnfold](https://hippunfold.readthedocs.io/en/latest/contributing/contributing.html).
-   GPU not required
-   Note: Apple M1 is currently **not supported**. We don't have a Docker arm64 container yet, and hippunfold is unusably slow with the emulated amd64 container. 

### Notes:

-   Inputs to Hippunfold should typically be a BIDS dataset including T1w images or T2w images. Higher-resolution data are preferred (\<= 0.8mm) but the pipeline will still work with 1mm T1w images. See [Tutorials](https://hippunfold.readthedocs.io/en/latest/tutorials/standardBIDS.html).
-   Other 3D imaging modalities (eg. ex-vivo MRI, 3D histology, etc.) can be used, but may require manual tissue segmentation as the current workflow relies on U-net segmentation trained only on common MRI modalities.


## Overview

There are several different ways of running HippUnfold. In order of increasing complexity/flexibility, we have:

1. CBRAIN Web-based Platform
2. Singularity Container on Linux
3. Docker Container on Windows/Mac (Intel)/Linux
4. Python Environment with Singularity Dependencies

### CBRAIN Web-based Platform

HippUnfold is available on the [CBRAIN plaform](https://github.com/aces/cbrain/wiki), a 
web-based platform for batch high-performance computing that is free for researchers.

Pros:
 - No software installation required
 - Fully point and click interface (no CLI)
 - Can perform batch-processing
Cons:
 - Must upload data for processing
 - Limited command-line options exposed
 - Cannot edit code


### Docker on Windows/Mac (Intel)/Linux

The HippUnfold BIDS App is available on a DockerHub as versioned releases and development branches.

Pros:
 - Compatible with non-Linux systems 
 - All dependencies+models in a single container
Cons:
 - Typically not possible on shared machines
 - Cannot use Snakemake cluster execution profiles
 - Cannot edit code

### Singularity Container

The same docker container can also be used with Singularity (now Apptainer). Instructions can be found below.

Pros:
 - All dependencies+models in a single container
 - Container stored as a single file (.sif)
Cons:
 - Compatible on shared systems with Singularity installed
 - Cannot use Snakemake cluster execution profiles
 - Cannot edit code


### Python Environment with Singularity Dependencies

Instructions for this can be found in the **Contributing** documentation page.

Pros:
 - Complete flexibility to modify code
 - External (Non-Python) dependencies as Singularity containers
Cons:
 - Must use Python virtual environment
 - Only compatible on Linux systems with Singularity for external dependencies



## Running HippUnfold with Docker


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
    /bids /output participant -np --modality T1w

Run it with maximum number of cores:

    docker run -it --rm \
    -v PATH_TO_BIDS_DIR:/bids:ro \
    -v PATH_TO_OUTPUT_DIR:/output \
    khanlab/hippunfold:latest \
    /bids /output participant -p --modality T1w --cores all

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

## Running HippUnfold with Singularity

Pull from dockerhub:

    singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest

See HippUnfold usage docs:

    singularity run -e khanlab_hippunfold_latest.sif -h

Do a dry run, printing the command at each step:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np --modality T1w

Run it with maximum number of cores:

    singularity run -e khanlab_hippunfold_latest.sif \
    PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -p --cores all --modality T1w

Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) to ensure this container can read and write to yout input and output directories, repsectively. For example, if your home directory is full or inaccessible, you may wish to set the following singularity parameters:

    export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity
    export SINGULARITY_BINDPATH=/YOURDIR:/YOURDIR

, where `YOURDIR` is your preferred storage location.
