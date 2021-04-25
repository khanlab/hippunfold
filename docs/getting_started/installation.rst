Installation
============

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and subfield segmentation)

Requirements
------------

* Docker (Mac/Windows/Linux) or Singularity (Linux)

* For the T2w (default) workflow: BIDS dataset with T1w and T2w images. 3D T2w TSE images with 0.8mm isotropic or higher resolution images are highly recommended. The T1w images are only used for linear registration, with the T2w images used for the convolutional neural network segmentation.

* GPU not required


Notes:
^^^^^^

#. T1w-only workflow is available too with ``--modality T1w``\ , however, this is discouraged unless you have high resolution (~0.7mm or better) T1w data, and performance will still likely be sub-optimal. This is currently being evaluated more thoroughly.

#. dwi workflows are also available but currently experimental



Running with Docker
-------------------

Pull the container::

   docker pull khanlab/hippunfold:latest

do a dry run, printing the command at each step::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -np 

run it with maximum number of cores::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -p --cores all


Running with Singularity
------------------------


If you are using singularity, there are two different ways to run ``hippunfold``. 

1. Pull the fully-contained BIDS app container from ``docker://khanlab/hippunfold:latest`` and run it

2. Clone the github repo, pip install, and run ``hippunfold``, which will pull the required singularity containers as needed.

If you want to modify the workflow at all or run a development branch, then option 2 is required. 

Furthermore, if you want to run the workflow on a Cloud system that Snakemake supports, or with a Cluster Execution profile, you should use Option 2 (see https://github.com/snakemake-profiles/doc and https://snakemake.readthedocs.io/en/stable/executing/cloud.html for more details on cluster profiles and cloud execution respectively).


Option 1:
^^^^^^^^^

Pull the container::
   
   singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest

do a dry run, printing the command at each step::

   singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 

run it with maximum number of cores::

   singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -p --cores all


Option 2:
^^^^^^^^^

This assumes you have created and activated a virtualenv or conda environment first.

#. Clone the repository::
   
    git clone --recursive http://github.com/khanlab/hippunfold

#. Install hippunfold using pip, with the -e option (for development mode)::

    pip install -e ./hippunfold

#. Run the following to download the U-net models::

    hippunfold_download_models

#. Now hippunfold will be installed for you and can perform a dry-run with::

    hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np

#. To run on all cores, and have snakemake pull any required containers, use::
    
    hippunfold  PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant --cores all --use-singularity




Additional instructions for Compute Canada 
------------------------------------------

Setting up a dev environment on graham:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are some instructions to get your python environment set-up on graham to run hippunfold:

#. Create a virtualenv and activate it::

      mkdir $SCRATCH/hippdev
      cd $SCRATCH/hippdev
      module load python/3
      virtualenv venv
      source venv/bin/activate

#. Follow the steps above to install from github repository

Running hippunfold jobs on graham:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In an interactive job (for testing)::
    
    regularInteractive -n 8
    hippunfold bids_dir out_dir participant --participant_label CC110037 -j 8


Submitting a job (for larger cores, more subjects), still single job, but snakemake will parallelize over the 32 cores::

    regularSubmit -j Fat hippunfold bids_dir out_dir participant  -j 32


Scaling up to ~hundred subjects (needs cc-slurm snakemake profile installed), submits 1 16core job per subject::
    
    hippunfold bids_dir out_dir participant  --profile cc-slurm


Scaling up to even more subjects (uses group-components to bundle multiple subjects in each job), 1 32core job for N subjects (e.g. 10)::
    
    hippunfold bids_dir out_dir participant  --profile cc-slurm --group-components subj=10

Note that this requires `neuroglia-helpers <https://github.com/khanlab/neuroglia-helpers>`_ for regularSubmit or regularInteractive wrappers, and the `cc-slurm <https://github.com/khanlab/cc-slurm>`_ snakemake profile for graham cluster execution with slurm. 
