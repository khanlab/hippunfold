# Contributing to Hippunfold

Hippunfold python package dependencies are managed with Poetry, which you\'ll need
installed on your machine. You can find instructions on the [poetry
website](https://python-poetry.org/docs/master/#installation).

HippUnfold also has a number of dependencies outside of python, including popular neuroimaging tools like `wb_command`, `ANTs`, `c3d`, and others listed in https://github.com/khanlab/hippunfold_deps. Thus we strongly recommend running HippUnfold with the `--use-singularity` flag, which will pull this container automatically and use it when required, unless you are comfortable using all of these tools yourself. 

Note: These instructions are only recommended if you are making changes to HippUnfold code to commit back to this repository, or if you are using Snakemake's cluster execution profiles. If not, it is easier to run HippUnfold when it is packaged into a single singularity container (e.g. `docker://khanlab/hippunfold:latest`). 


## Set-up your development environment:

Clone the repository and install dependencies and dev dependencies with
poetry:

    git clone http://github.com/khanlab/hippunfold
    cd hippunfold
    poetry install

Poetry will automatically create a virtual environment. To customize where 
these virtual environments are stored see poetry docs 
[here](https://python-poetry.org/docs/configuration/)

Then, you can run hippunfold with:

    poetry run hippunfold

or you can activate a virtualenv shell and then run hippunfold directly:

    poetry shell
    hippunfold

You can exit the poetry shell with `exit`.

Note: you can alternatively use `pip install` to install dependencies.

## Running code format quality checking and fixing:

Hippunfold uses [poethepoet](https://github.com/nat-n/poethepoet) as a
task runner. You can see what commands are available by running:

    poetry run poe

We use `black` and `snakefmt` to ensure
formatting and style of python and Snakefiles is consistent. There are
two task runners you can use to check and fix your code, and can be
invoked with:

    poetry run poe quality_check
    poetry run poe quality_fix

Note that if you are in a poetry shell, you do not need to prepend
`poetry run` to the command.

## Dry-run testing your workflow:

Using Snakemake\'s dry-run option (`--dry-run`/`-n`) is an easy way to verify any
changes to the workflow are working correctly. The `test_data` folder contains a 
number of *fake* bids datasets (i.e. datasets with zero-sized files) that are useful
for verifying different aspects of the workflow. These dry-run tests are
part of the automated github actions that run for every commit.

You can use the hippunfold CLI to perform a dry-run of the workflow,
e.g. here printing out every command as well:

    hippunfold test_data/bids_singleT2w test_out participant --modality T2w --use-singularity -np

As a shortcut, you can also use `snakemake` instead of the
hippunfold CLI, as the `snakebids.yml` config file is set-up
by default to use this same test dataset, as long as you run snakemake
from the `hippunfold` folder that contains the
`workflow` folder:

    cd hippunfold
    snakemake -np

## Wet-run testing your workflow:

Opensource data has been made available for wet-run testing any changes to HippUnfold on OSF [here](https://osf.io/k2nme/). These are real opensource data meant ot span a wide array of HippUnfold use cases with various modalities, resolutions, developmental stages, and species. Please note that some changes to hippUnfold will impact all worflows and should be run and then visually inspected for all these cases. Other changes may impact only one workflow (e.g. adding a new template) and so the recommended run parameters should be adjusted accordingly. 

## Instructions for Compute Canada

This section provides an example of how to set up a `pip installed` copy
of HippUnfold on Compute Canada\'s `graham` cluster.

### Setting up a dev environment on graham:

Here are some instructions to get your python environment set-up on
graham to run HippUnfold:

1.  Create a virtualenv and activate it:

        mkdir $SCRATCH/hippdev
        cd $SCRATCH/hippdev
        module load python/3.8
        virtualenv venv
        source venv/bin/activate

2.  Install HippUnfold

        git clone https://github.com/khanlab/hippunfold.git
        pip install hippunfold/
        
3. To run Hippunfold on Graham as a member of the Khan lab, please configure the
[neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) with the khanlab profile.

4. To avoid having to download trained models (see section [below](#deep-learning-nnu-net-model-files)), you can set an environment variable in your bash profile (~/.bash_profile) with the location of the
trained models. For Khan lab's members, the following line must be add to the bash profile file:

        export HIPPUNFOLD_CACHE_DIR="/project/6050199/akhanf/opt/hippunfold_trained_models"


Note: make sure to reload your bash profile if needed (`source ~./bash_profile`).        

5. For an easier execution in Graham, it's recommended to also install
[cc-slurm](https://github.com/khanlab/cc-slurm) snakemake profile for cluster execution with slurm.

Note if you want to run hippunfold with modifications to your cloned 
repository, you either need to pip install again, or run hippunfold the following, since 
an `editable` pip install is not allowed with pyproject:

        python <YOUR_HIPPUNFOLD_DIR>/hippunfold/run.py

### Running hippunfold jobs on graham:

Note that this requires
[neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) for
regularSubmit or regularInteractive wrappers, and the
[cc-slurm](https://github.com/khanlab/cc-slurm) snakemake profile for cluster execution with slurm.

In an interactive job (for testing):

    regularInteractive -n 8
    hippunfold <PATH_TO_BIDS_DIR> <PATH_TO_OUTPUT_DIR> participant \
    --participant_label 001 -j 8 --modality T1w --use-singularity \
    --singularity-prefix $SNAKEMAKE_SINGULARITY_DIR

Where:
 - `--participant_label 001` is used to specify only one subject from a BIDS
directory presumeably containing many subjects.
 - `-j 8` specifies the number of cores used
 - `--modality T1w` is used to specify that a T1w dataset is being processed
 - `--singularity-prefix $SNAKEMAKE_SINGULARITY_DIR` specifies the directory in
which singularity images will be stored. The environment variable is created
when installing neuroglia-helpers.

Submitting a job (for larger cores, more subjects), still single job,
but snakemake will parallelize over the 32 cores:

    regularSubmit -j Fat \
    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant  -j 32 \
    --modality T1w --use-singularity --singularity-prefix $SNAKEMAKE_SINGULARITY_DIR

Scaling up to \~hundred subjects (needs cc-slurm snakemake profile
installed), submits 1 16core job per subject:

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant \
    --modality T1w --use-singularity --singularity-prefix $SNAKEMAKE_SINGULARITY_DIR \
    --profile cc-slurm

Scaling up to even more subjects (uses group-components to bundle
multiple subjects in each job), 1 32core job for N subjects (e.g. 10):

    hippunfold PATH_TO_BIDS_DIR PATH_TO_OUTPUT_DIR participant \
    --modality T1w --use-singularity --singularity-prefix $SNAKEMAKE_SINGULARITY_DIR \
    --profile cc-slurm --group-components subj=10

### Running hippunfold jobs on the CBS server
1. Clone the repository and install dependencies and dev dependencies with poetry:

       git clone http://github.com/khanlab/hippunfold
       cd hippunfold
       poetry install
If poetry is not installed, please refer to the [installation documentation](https://python-poetry.org/docs/). If the command poetry is not found, add the following line to your bashrc file located in your home directory (considering that the poetry binary is located under `$HOME/.local/bin`:

       export PATH=$PATH:$HOME/.local/bin
2. To avoid having to download containers and trained models (see section [below](#deep-learning-nnu-net-model-files)), add the `$SNAKEMAKE_SINGULARITY_DIR` and `$HIPPUNFOLD_CACHE_DIR` environment variables to the bashrc file. For Khan lab's members, add the following lines:

        export SNAKEMAKE_SINGULARITY_DIR="/cifs/khan/shared/containers/snakemake_containers"
        export HIPPUNFOLD_CACHE_DIR="/cifs/khan/shared/data/hippunfold_models"

3. HippUnfold might be executed using `poetry run hippunfold <arguments>` or through the `poetry shell` method. Refer to previous section for more information in regards to execution options. 

4. On the CBS server you should always set your output folder to a path inside `/localscratch`, and not your home folder or a `/srv` or `/cifs` path, and copy the final results out after they have finished computing. Please be aware that the CBS server may not be the most efficient option for running a large number of subjects (since you are limited in processing cores vs a HPC cluster).

5. If you are using input files in your home directory (or in your `graham` mount in your home directory), you may also need to also add the following to your bashrc file (Note: this will become a default system-enabled option soon)

        export SINGULARITY_BINDPATH="/home/ROBARTS:/home/ROBARTS"

## Deep learning nnU-net model files

The trained model files we use for hippunfold are large and thus are not
included directly in this github repository, and instead are downloaded
from Zenodo releases. 

### For HippUnfold versions earlier than 1.3.0 (< 1.3.0): 
If you are using the docker/singularity container, `docker://khanlab/hippunfold`, they are pre-downloaded there, in `/opt/hippunfold_cache`.

If you are not using this container, you will need to download the models before running hippunfold, by running:

    hippunfold_download_models
    
This console script (installed when you install hippunfold) downloads all the models to a cache dir on your system, 
which on Linux is typically `~/.cache/hippunfold`. To override this, you can set the `HIPPUNFOLD_CACHE_DIR` environment
variable before running `hippunfold_download_models` and `hippunfold`.

### NEW: For HippUnfold versions 1.3.0 and later (>= 1.3.0):
With the addition of new models, including all models in the container was not feasible and a change was made to 
**not include** any models in the docker/singularity containers. In these versions, the `hippunfold_download_models` command
is removed, and any models will simply be downloaded as part of the workflow. As before, all models will be stored in the system cache dir, 
which is typically `~/.cache/hippunfold`, and to override this can set the `HIPPUNFOLD_CACHE_DIR` environment variable before running `hippunfold`.

If you want to pre-download a model (e.g. if your compute nodes do not have internet access), you can run simply run `download_model` rule in HippUnfold e.g.:

```
hippunfold BIDS_DIR OUTPUT_DIR PARTICIPANT_LEVEL --modality T1w --until download_model -c 1
```


## Overriding Singularity cache directories

By default, singularity stores image caches in your home directory when you run `singularity pull` or `singularity run`. As described above, hippunfold also stores deep learning models in your home directory. If your home directory is full or otherwise inaccessible, you may want to change this with the following commands:

    export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity
    export SINGULARITY_BINDPATH=/YOURDIR:/YOURDIR
    export HIPPUNFOLD_CACHE_DIR=/YOURDIR/.cache/hippunfold/
    
If you are running `hippunfold` with the `--use-singularity` option, hippunfold will download the required singularity containers for rules that require it. These containers are placed in the `.snakemake` folder in your hippunfold output directory, but this can be overriden with the Snakemake option: `--singularity-prefix DIRECTORY`
  

