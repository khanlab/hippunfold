Hippunfold
==========

.. image:: https://images.microbadger.com/badges/version/khanlab/hippunfold.svg
   :target: https://microbadger.com/images/khanlab/hippunfold "Get your own version badge on microbadger.com"
   :alt: 


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


