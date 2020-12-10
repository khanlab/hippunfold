# hippunfold
BIDS App for Hippocampal AutoTop (automated hippocampal unfolding and subfield segmentation)

## requires:

Singularity
Python (>=3.7)

## install:

0. Clone this repo
```
git clone http://github.com/khanlab/hippunfold
```

1. Get dependent containers and update paths in config/snakebids.yml (TODO: incoporate autodownload into snakebids, pip)
```
docker://khanlab/prepdwi:v0.0.13
docker://khanlab/autotop_deps:v0.1
docker://kaczmarj/ants:2.2.0
```

(plus FSL cuda container required if running dwi preprocessing)

2. Pip install

`pip install -e <path to cloned hippunfold>`



## setting up dev environemnt
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
git clone http://github.com/akhanf/snakebids
git clone --recursive http://github.com/khanlab/hippunfold
```

3. install snakebids using pip, with the -e option (for development mode):
```
pip install -e ./snakebids
pip install -e ./hippunfold
```

Now hippunfold will be installed for you and can run with:

hippunfold  <args here> 

Any containers used are included in the hippunfold workflow, and if in khanlab group on graham, will already be good to go..  If you log out, you just need to re-activate the virtualenv to start again 

If you ever want the latest code, can just pull it:
```
cd snakebids
git pull
```

or if you need a branch, can: `git checkout <name of branch>`
```


## running on graham:

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
