[![Documentation Status](https://readthedocs.org/projects/hippunfold/badge/?version=latest)](https://hippunfold.readthedocs.io/en/latest/?badge=latest)
[![CircleCI](https://circleci.com/gh/khanlab/hippunfold.svg?style=svg)](https://circleci.com/gh/khanlab/hippunfold)
![PyPI](https://img.shields.io/pypi/v/hippunfold)
![Docker Pulls](https://img.shields.io/docker/pulls/khanlab/hippunfold)

# Hippunfold

This tool aims to automatically model the topological folding structure
of the human hippocampus, and computationally unfold the hippocampus to 
segment subfields and generate hippocampal and dentate gyrus surfaces.

The overall workflow can be summarized in the following steps:

![Pipeline Overview](https://raw.githubusercontent.com/khanlab/hippunfold/master/docs/images/hippunfold_overview.jpg)

1.  Pre-processing MRI images including non-uniformity correction, 
    resampling to 0.3mm isotropic subvolumes, registration and cropping to coronal-oblique 
    subvolumes around each hippocampus
2.  Automatic segmentation of hippocampal tissues and surrounding
    structures via deep convolutional neural network U-net (nnU-net), models available
    for T1w, T2w, hippocampal b500 dwi, and neonatal T1w, and post-processing
    with fluid label-label registration to a high-resolution average template
3.  Imposing of coordinates across the anterior-posterior and 
    proximal-distal dimensions of hippocampal grey matter
    via solving Laplace's equation, and using an equivolume solution for
    laminar coordinates
4.  Generating transformations between native and unfolded spaces using scattered 
    interpolation on the hippocampus and dentate gyrus coordinates separately
5.  Applying these transformations to generate surfaces meshes of the hippocampus
    and dentate gyrus, and extraction of morphological surface-based features including
    thickness, curvature, and gyrification index, sampled at the midthickness surface, and mapping 
    subfield labels from the histological BigBrain atlas of the hippocampus 
6.  Generating high-resolution volumetric segmentations of the subfields using the same
    transformations and volumetric representations of the coordinates.


**Full Documentation:**  [here](https://hippunfold.readthedocs.io**)


**Relevant papers:**
-  DeKraker J, Haast RAM, Yousif MD, Karat B, Köhler S, Khan AR. HippUnfold: Automated hippocampal unfolding, morphometry, and subfield segmentation. bioRxiv 2021.12.03.471134; doi: 10.1101/2021.12.03.471134 [link](https://www.biorxiv.org/content/10.1101/2021.12.03.471134v1)
- DeKraker J, Ferko KM, Lau JC, Köhler S, Khan AR. Unfolding the hippocampus: An intrinsic coordinate system for subfield segmentations and quantitative mapping. Neuroimage. 2018 Feb 15;167:408-418. doi: 10.1016/j.neuroimage.2017.11.054. Epub 2017 Nov 23. PMID: 29175494. [link](https://pubmed.ncbi.nlm.nih.gov/29175494/)
- DeKraker J, Lau JC, Ferko KM, Khan AR, Köhler S. Hippocampal subfields revealed through unfolding and unsupervised clustering of laminar and morphological features in 3D BigBrain. Neuroimage. 2020 Feb 1;206:116328. doi: 10.1016/j.neuroimage.2019.116328. Epub 2019 Nov 1. PMID: 31682982. [link](https://pubmed.ncbi.nlm.nih.gov/31682982/)
- DeKraker J, Köhler S, Khan AR. Surface-based hippocampal subfield segmentation. Trends Neurosci. 2021 Nov;44(11):856-863. doi: 10.1016/j.tins.2021.06.005. Epub 2021 Jul 22. PMID: 34304910. [link](https://pubmed.ncbi.nlm.nih.gov/34304910/)


