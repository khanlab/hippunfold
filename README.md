[![Documentation Status](https://readthedocs.org/projects/hippunfold/badge/?version=latest)](https://hippunfold.readthedocs.io/en/latest/?badge=latest)
![Docker Pulls](https://img.shields.io/docker/pulls/khanlab/hippunfold)
![Version](https://img.shields.io/github/v/tag/khanlab/hippunfold?label=version)

**Full Documentation:**  [here](https://hippunfold.readthedocs.io/en/latest/?badge=latest)

# Hippunfold

This tool aims to automatically model the topological folding structure
of the human hippocampus, and computationally unfold the hippocampus to 
segment subfields and generate hippocampal and dentate gyrus surfaces.

![Hippo Fold Unfold](images/subfields_foldunfold.png)

This is especially useful for:
- Visualization
- Topologically-constrained intersubject registration
- Parcellation (ie. registration to an unfolded atlas)
- Morphometry (eg. thickness, surface area, curvature, and gyrification measures)
- Quantitative mapping (eg. map your qT1 MRI data to a midthickness surface; extract laminar profiles perpendicular to this surface)

## NEW: Version 1.3.0 release 

Major changes include the addition of unfolded space registration to a reference atlas harmonized across 7 ground-truth histology samples. This method allows shifting in unfolded space, providing even better intersubject alignment.

Read more in our [preprinted manuscript](https://www.biorxiv.org/content/10.1101/2023.03.30.534978v1)

Also the ability to specify a new **experimental** UNet model that is contrast-agnostic using [synthseg](https://github.com/BBillot/SynthSeg) and trained using more detailed segmentations. This generally produces more detailed results but has not been extensively tested yet. 

## Workflow

The overall workflow can be summarized in the following steps:

![Pipeline Overview](images/hippunfold_overview_unfoldreg.png)

For more information, see
**Full Documentation:**  [here](https://hippunfold.readthedocs.io/en/latest/?badge=latest)

## Additional tools

For plotting, mapping fMRI, DWI or other data, and manipulating surfaces, see [here](https://github.com/jordandekraker/hippunfold_toolbox)

For statistical testing (spin tests) in unfolded space, see [here](https://github.com/Bradley-Karat/Hippo_Spin_Testing)

## Relevant papers:
- DeKraker, J., Haast, R. A., Yousif, M. D., Karat, B., Lau, J. C., Köhler, S., & Khan, A. R. (2022). Automated hippocampal unfolding for morphometry and subfield segmentation with HippUnfold. Elife, 11, e77945. [link](https://doi.org/10.7554/eLife.77945)
- DeKraker, J., Palomero-Gallagher, N., Kedo, O., Ladbon-Bernasconi, N., Muenzing, S. E., Axer, M., ... & Evans, A. C. (2023). Evaluation of surface-based hippocampal registration using ground-truth subfield definitions. bioRxiv, 2023-03. [link](https://www.biorxiv.org/content/10.1101/2023.03.30.534978v1)
- DeKraker J, Ferko KM, Lau JC, Köhler S, Khan AR. Unfolding the hippocampus: An intrinsic coordinate system for subfield segmentations and quantitative mapping. Neuroimage. 2018 Feb 15;167:408-418. doi: 10.1016/j.neuroimage.2017.11.054. Epub 2017 Nov 23. PMID: 29175494. [link](https://pubmed.ncbi.nlm.nih.gov/29175494/)
- DeKraker J, Lau JC, Ferko KM, Khan AR, Köhler S. Hippocampal subfields revealed through unfolding and unsupervised clustering of laminar and morphological features in 3D BigBrain. Neuroimage. 2020 Feb 1;206:116328. doi: 10.1016/j.neuroimage.2019.116328. Epub 2019 Nov 1. PMID: 31682982. [link](https://pubmed.ncbi.nlm.nih.gov/31682982/)
- DeKraker J, Köhler S, Khan AR. Surface-based hippocampal subfield segmentation. Trends Neurosci. 2021 Nov;44(11):856-863. doi: 10.1016/j.tins.2021.06.005. Epub 2021 Jul 22. PMID: 34304910. [link](https://pubmed.ncbi.nlm.nih.gov/34304910/)


