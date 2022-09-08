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

HippUnfold is available on the [CBRAIN platform](https://github.com/aces/cbrain/wiki), a 
web-based platform for batch high-performance computing that is free for researchers.

#### Pros:
- No software installation required
- Fully point and click interface (no CLI)
- Can perform batch-processing

#### Cons:
- Must upload data for processing
- Limited command-line options exposed
- Cannot edit code


### Docker on Windows/Mac (Intel)/Linux

The HippUnfold BIDS App is available on a DockerHub as versioned releases and development branches.

#### Pros:
- Compatible with non-Linux systems 
- All dependencies+models in a single container

#### Cons:
- Typically not possible on shared machines
- Cannot use Snakemake cluster execution profiles
- Cannot edit code

### Singularity Container

The same docker container can also be used with Singularity (now Apptainer). Instructions can be found below.

#### Pros:
- All dependencies+models in a single container
- Container stored as a single file (.sif)

#### Cons:
- Compatible on shared systems with Singularity installed
- Cannot use Snakemake cluster execution profiles
- Cannot edit code


### Python Environment with Singularity Dependencies

Instructions for this can be found in the **Contributing** documentation page.

#### Pros:
- Complete flexibility to modify code
- External (Non-Python) dependencies as Singularity containers

#### Cons:
- Must use Python virtual environment
- Only compatible on Linux systems with Singularity for external dependencies



## Running HippUnfold with Docker


Note: These instructions assume you have Docker installed already on your system.

Download and extract a single-subject BIDS dataset for this test:

    wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
    tar -xvf hippunfold_test_data.tar

This will create a `ds002168/` folder with a single subject, that has a 
both T1w and T2w images. 


Pull the container:

    docker pull khanlab/hippunfold:latest

Run HippUnfold without any arguments to print the short help:

    docker run -it --rm \
    khanlab/hippunfold:latest    

Use the `-h` option to get a detailed help listing:

    docker run -it --rm \
    khanlab/hippunfold:latest \
    -h

Note that all the Snakemake command-line options are also available in
HippUnfold, and can be listed with `--help-snakemake`:

    docker run -it --rm \
    khanlab/hippunfold:latest \
    --help-snakemake


Now let's run it on the test dataset. The `--modality` flag is a 
required argument, and describes what image we use for segmentation. Here 
we will use the T1w image. We will also use the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -n


For those not familiar with Docker, the first three lines of this
example are generic Docker arguments to ensure it is run with the safest
options and has permission to access your input and output directories
The fourth line specifies the HippUnfold Docker container, and the fifth line contains the required
arguments for HippUnfold, after which you can additionally specify optional arguments. You may want to familiarize yourself with
[Docker options](https://docs.docker.com/engine/reference/run/), and an
overview of HippUnfold arguments is provided in the [Command line
interface](https://hippunfold.readthedocs.io/en/latest/usage/cli.html)
documentation section.


The first three arguments to HippUnfold (as with any BIDS App) are the input
folder, the output folder, and then the analysis level. The `participant` analysis 
level is used in HippUnfold for performing the segmentation, unfolding, and any
participant-level processing. The `group` analysis is used to combine subfield volumes
across subjects into a single tsv file.


When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -np | less


Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells HippUnfold how many cores to use.
 Using `--cores 8` means that HippUnfold will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (hippunfold on a single subject) may take ~30 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -p --cores all


After this completes, you should have a `ds002168_hippunfold` folder with outputs for the one subject.

If you alternatively want to run HippUnfold using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. In this case, since the T2w image in the 
test dataset has a limited FOV, we should also make use of the `--t1-reg-template` command-line option,
which will make use of the T1w image for template registration, since a limited FOV T2w template does not exist.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold_t2w:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T2w --t1-reg-template -p --cores all

Note that if you run with a different modality, you should use a separate output folder, since some of the files 
would be overwritten if not.


## Running HippUnfold with Singularity

Note: These instructions assume you have Singularity installed on your system.

Download and extract a single-subject BIDS dataset for this test:

    wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
    tar -xvf hippunfold_test_data.tar

This will create a `ds002168/` folder with a single subject, that has a 
both T1w and T2w images. 


Next we will pull the container from dockerhub:

    singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest


If you encounter any errors  pulling the container, it may be because you are running 
out of disk space in your cache folders. Note, you can change these locations 
by setting environment variables, e.g.:
    
    export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity


Run HippUnfold any arguments to print the short help:

    singularity run -e khanlab_hippunfold_latest.sif 

Use the `-h` option to get a detailed help listing:

    singularity run -e khanlab_hippunfold_latest.sif -h

Note that all the Snakemake command-line options are also available in
HippUnfold, and can be listed with `--help-snakemake`:

    singularity run -e khanlab_hippunfold_latest.sif --help-snakemake


Now let's run it on the test dataset. The `--modality` flag is a 
required argument, and describes what image we use for segmentation. Here 
we will use the T1w image. We will also use the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.


    singularity run -e khanlab_hippunfold_latest.sif \
    ds002168 ds002168_hippunfold participant -n --modality T1w


The first three arguments to HippUnfold (as with any BIDS App) are the input
folder, the output folder, and then the analysis level. The `participant` analysis 
level is used in HippUnfold for performing the segmentation, unfolding, and any
participant-level processing. The `group` analysis is used to combine subfield volumes
across subjects into a single tsv file.


When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    singularity run -e khanlab_hippunfold_latest.sif \
    ds002168 ds002168_hippunfold participant -np --modality T1w | less


Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells HippUnfold how many cores to use.
 Using `--cores 8` means that HippUnfold will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (hippunfold on a single subject) may take ~30 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.


    singularity run -e khanlab_hippunfold_latest.sif \
    ds002168 ds002168_hippunfold participant -p --cores all --modality T1w


Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) to ensure the container can read and write to yout input and output directories, respectively. You can bind paths easily by setting an 
environment variable, e.g. if you have a `/project` folder that contains your data, you can add it to the `SINGULARITY_BINDPATH` so it is available when you are running a container:

    export SINGULARITY_BINDPATH=/data:/data



After this completes, you should have a `ds002168_hippunfold` folder with outputs for the one subject.

If you alternatively want to run HippUnfold using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. In this case, since the T2w image in the 
test dataset has a limited FOV, we should also make use of the `--t1-reg-template` command-line option,
which will make use of the T1w image for template registration, since a limited FOV T2w template does not exist.

    singularity run -e khanlab_hippunfold_latest.sif \
    ds002168 ds002168_hippunfold_t2w participant --modality T2w --t1-reg-template -p --cores all

Note that if you run with a different modality, you should use a separate output folder, since some of the files 
would be overwritten if not.




