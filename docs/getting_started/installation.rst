Installation
============

BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and subfield segmentation)

Requirements
------------

Docker (Mac/Windows/Linux) or Singularity (Linux)

BIDS dataset with T1w and T2w images. Highly recommend using 0.8mm isotropic or higher 3D T2w TSE. 

Note 1: T1w only workflow is available too with ``--modality T1w``\ , however, this is discouraged unless you have high resolution (~0.7mm or better) T1w data, and performance will likely be sub-optimal. This is currently being evaluated.

Note 2: dwi workflows are also available but currently experimental


Docker:
^^^^^^^

Pull the container:

.. code-block::

   docker pull khanlab/hippunfold:latest

do a dry run, printing the command at each step:

.. code-block::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -np 

run it with maximum number of cores:

.. code-block::

   docker run -it --rm -v PATH_TO_BIDS_DIR:/bids:ro -v PATH_TO_OUTPUT_DIR:/output khanlab/hippunfold:latest /bids /output participant -p --cores all


Singularity:
^^^^^^^^^^^^

Pull the container:

.. code-block::

   singularity pull khanlab_hippunfold_latest.sif docker://khanlab/hippunfold:latest

do a dry run, printing the command at each step:

.. code-block::

   singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant -np 

run it with maximum number of cores:

.. code-block::

   singularity run -e khanlab_hippunfold_latest.sif khanlab/hippunfold:latest PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -p --cores all


Setting up a dev environment:
-----------------------------

Here are some instructions to get your python environment set-up on graham to run hippunfold:


#. create a virtualenv and activate it:
   .. code-block::

      mkdir $SCRATCH/hippdev
      cd $SCRATCH/hippdev
      module load python/3
      virtualenv venv
      source venv/bin/activate

#. 
   clone the source repos (so you can make/pull changes easily, or change to branch):

   .. code-block::

      git clone --recursive http://github.com/khanlab/hippunfold

#. 
   install snakebids using pip, with the -e option (for development mode):

   .. code-block::

      pip install -e ./hippunfold

Now hippunfold will be installed for you and can run with:

.. code-block::

   hippunfold  <args here> 


Any containers used are included in the hippunfold workflow, and if in khanlab group on graham, will already be good to go..  If you log out, you just need to re-activate the virtualenv to start again 

If you ever want the latest code, can just pull it:

.. code-block::

   cd hippunfold
   git pull

or if you need a branch, can: ``git checkout <name of branch>``

running on graham:
^^^^^^^^^^^^^^^^^^

In an interactive job (for testing):

.. code-block::

   regularInteractive -n 8
   hippunfold bids_dir out_dir participant --participant_label CC110037 -j 8


Submitting a job (for larger cores, more subjects), still single job, but snakemake will parallelize over the 32 cores:

.. code-block::

   regularSubmit -j Fat hippunfold bids_dir out_dir participant  -j 32


Scaling up to ~hundred subjects (needs cc-slurm snakemake profile installed), submits 1 16core job per subject:

.. code-block::

   hippunfold bids_dir out_dir participant  --profile cc-slurm


Scaling up to even more subjects (uses group-components to bundle multiple subjects in each job), 1 32core job for N subjects (e.g. 10):

.. code-block::

   hippunfold bids_dir out_dir participant  --profile cc-slurm --group-components subj=10



