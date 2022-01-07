Hippunfold
==========

.. image:: https://readthedocs.org/projects/hippunfold/badge/?version=latest
   :target: https://hippunfold.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. image:: https://circleci.com/gh/khanlab/hippunfold.svg?style=svg
   :target: https://circleci.com/gh/khanlab/hippunfold
   :alt: CircleCI



This tool aims to automatically model the topological folding structure of the human hippocampus. It is currently set up to use sub-millimetric T2w MRI data, but may be adapted for other data types. This can then be used to apply the hippocampal unfolding methods presented in `DeKraker et al., 2019 <https://www.sciencedirect.com/science/article/pii/S1053811917309977>`_, and ex-vivo subfield boundaries can be topologically applied from `DeKraker et al., 2020 <https://www.sciencedirect.com/science/article/pii/S105381191930919X?via%3Dihub>`_.

.. image:: https://github.com/khanlab/hippunfold/raw/master/docs/pipeline_overview.png
    :align: center
    :alt: Pipeline Overview

The overall workflow can be summarized in the following steps:

0. Resampling to a 0.3mm isotropic, coronal oblique, cropped hippocampal block

1. Automatic segmentation of hippocampal tissues and surrounding structures via deep convolutional neural network U-net `Li _et al_., 2017 <https://arxiv.org/abs/1707.01992>`_ _OR_ Manual segmentation of hippocampal tissues and surrounding structures using `this <https://ars.els-cdn.com/content/image/1-s2.0-S1053811917309977-mmc1.pdf>`_ protocol

2. Post-processing via fluid label-label registration to a high resolution, topoligically correct averaged template

3. Imposing of coordinates across the anterior-posterior, proximal-distal, and laminar dimensions of hippocampal grey matter via solving the Laplace equation

4. Extraction of a grey matter mid-surface and morpholigical features (thickness, curvature, gyrification index, and, if available, quantitative MRI values sampled along the mid-surface for reduced partial-voluming)

5. Quality assurance via inspection of Laplace gradients, grey matter mid-surface, and flatmapped features

6. Application of subfield boundaries according to predifined topological coordinates



Contributing
============

Hippunfold dependencies are managed with Poetry, which you'll need installed on your machine. You can find instructions on the `poetry website <https://python-poetry.org/docs/master/#installation>`_. 

Set-up your development environment:
------------------------------------

Clone the repository and install dependencies and dev dependencies with poetry::

   git clone http://github.com/khanlab/hippunfold
   cd hippunfold
   poetry install

Poetry will automatically create a virtualenv. To customize ... (TODO: finish this part)

Then, you can run hippunfold with::

   poetry run hippunfold
   
or you can activate a virtualenv shell and then run hippunfold directly::

   poetry shell
   hippunfold
   
You can exit the poetry shell with `exit`.


Running code format quality checking and fixing:
------------------------------------------------

Hippunfold uses `poethepoet <https://github.com/nat-n/poethepoet>`_ as a task runner. You can see what commands are available by running::

    poetry run poe
    
We use `black` and `snakefmt` to ensure formatting and style of python and Snakefiles is consistent. There are two task runners you can use to check and fix your code, and can be invoked with::

   poetry run poe quality_check
   poetry run poe quality_fix

Note that if you are in a poetry shell, you do not need to prepend `poetry run` to the command. 

Dry-run testing your workflow:
------------------------------

Using Snakemake's dry-run option (`--dry-run`/`-n`) is an easy way to verify any changes to the workflow are working correctly. The `test_data` folder contains a number of *fake* bids datasets (i.e. datasets with zero-sized files) that are useful for verifying different aspects of the workflow. These dry-run tests are part of the automated github actions that run for every commit. 

You can use the hippunfold CLI to perform a dry-run of the workflow, e.g. here printing out every command as well::

   hippunfold test_data/bids_singleT2w test_out participant --modality T2w -np

As a shortcut, you can also use `snakemake` instead of the hippunfold CLI, as the `snakebids.yml` config file is set-up by default to use this same test dataset, as long as you run snakemake from the `hippunfold` folder that contains the `workflow` folder::

   cd hippunfold
   snakemake -np
   


   
