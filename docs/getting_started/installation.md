# Installation

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and
subfield segmentation)

## Requirements

HippUnfold can be run either with Conda, or with containers
-   Conda (Linux/Mac only, no containers needed)
-   **--OR--**  Docker (Intel Mac/Windows/Linux) or Apptainer (formerly known as Singularity) (Linux)
-   For those wishing to contribute or modify the code, see [Contributing to HippUnfold](https://hippunfold.readthedocs.io/en/latest/contributing/contributing.html).
-   GPU not required

### Notes:

-   Inputs to HippUnfold should typically be a BIDS dataset including T1w images or T2w images. Higher-resolution data are preferred (\<= 0.8mm) but the pipeline will still work with 1mm T1w images. See [Tutorials](https://hippunfold.readthedocs.io/en/latest/tutorials/standardBIDS.html).
-   Other 3D imaging modalities (eg. ex-vivo MRI, 3D histology, etc.) can be used, but may require manual tissue segmentation as the current workflow relies on U-net segmentation trained only on common MRI modalities.


## Comparison of methods for running HippUnfold

There are several different ways of running HippUnfold:

1. Conda Environment (Linux/macOS)
2. Apptainer/Singularity Container on Linux
3. Docker Container on Windows/Mac (Intel)/Linux
4. CBRAIN Web-based Platform

### Conda

As of version **2.0.0**, HippUnfold is available via [Conda](https://docs.conda.io/), offering a container-free way to run the tool on Linux and macOS systems.

#### Pros:
- No need for any containers
- Easy to install via `conda`
- Compatible with Snakemake execution profiles
- Good option for users who prefer Python virtual environments

#### Cons:
- Not compatible with Windows
- Requires installation of Conda or Miniconda

### Apptainer

The same docker container can also be used with Singularity (now Apptainer). Instructions can be found below.

#### Pros:
- All dependencies+models (* See Note 1) in a single container 
- Container stored as a single file (.sif)

#### Cons:
- Compatible on shared systems with Singularity installed
- Cannot use Snakemake cluster execution profiles
- Cannot edit code

### Docker

The HippUnfold BIDS App is available on a DockerHub as versioned releases and development branches.

#### Pros:
- Compatible with non-Linux systems 
- All dependencies+models (* See Note 1) in a single container

#### Cons:
- Typically not possible on shared machines
- Cannot use Snakemake cluster execution profiles
- Cannot edit code


### CBRAIN

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


## Note 1: 
As of version 1.3.0 of HippUnfold, containers are no longer shipped with all the models, and the models are downloaded as part of the workflow. By default, models are placed in `~/.cache/hippunfold` unless you set the `HIPPUNFOLD_CACHE_DIR` environment variable. See [Deep learning nnU-net model files](https://hippunfold.readthedocs.io/en/latest/contributing/contributing.html#deep-learning-nnu-net-model-files) for more information.

