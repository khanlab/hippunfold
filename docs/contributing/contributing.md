# Contributing to Hippunfold

Hippunfold dependencies are managed with Poetry, which you\'ll need
installed on your machine. You can find instructions on the [poetry
website](https://python-poetry.org/docs/master/#installation).

## Set-up your development environment:

Clone the repository and install dependencies and dev dependencies with
poetry:

    git clone http://github.com/khanlab/hippunfold
    cd hippunfold
    poetry install

Poetry will automatically create a virtualenv. To customize \... (TODO:
finish this part)

Then, you can run hippunfold with:

    poetry run hippunfold

or you can activate a virtualenv shell and then run hippunfold directly:

    poetry shell
    hippunfold

You can exit the poetry shell with `exit`.

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
changes to the workflow are working correctly. The `test_data`
 folder contains a number of *fake* bids
datasets (i.e. datasets with zero-sized files) that are useful for
verifying different aspects of the workflow. These dry-run tests are
part of the automated github actions that run for every commit.

You can use the hippunfold CLI to perform a dry-run of the workflow,
e.g. here printing out every command as well:

    hippunfold test_data/bids_singleT2w test_out participant --modality T2w -np

As a shortcut, you can also use `snakemake` instead of the
hippunfold CLI, as the `snakebids.yml` config file is set-up
by default to use this same test dataset, as long as you run snakemake
from the `hippunfold` folder that contains the
`workflow` folder:

    cd hippunfold
    snakemake -np
