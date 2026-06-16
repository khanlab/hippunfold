[![Documentation Status](https://readthedocs.org/projects/hippunfold/badge/?version=latest)](https://hippunfold.readthedocs.io/en/latest/?badge=latest)
![Docker Pulls](https://img.shields.io/docker/pulls/khanlab/hippunfold)
![Version](https://img.shields.io/github/v/tag/khanlab/hippunfold?label=version)
<img align="right" width="200" src="https://github.com/khanlab/hippunfold/assets/25106300/0c16d33e-893a-4ac3-b127-21fa843823d5">

**Full Documentation:**  [here](https://hippunfold.readthedocs.io/en/latest/?badge=latest)

# HippUnfold

This tool aims to automatically model the topological folding structure
of the human hippocampus, and computationally unfold it.
![Hippo Fold Unfold](https://raw.githubusercontent.com//khanlab/hippunfold/master/docs/images/subfields_foldunfold.png)

This is especially useful for:
- Visualization
- Topologically-constrained intersubject registration
- Parcellation (ie. registration to an unfolded atlas)
- Morphometry (eg. thickness, surface area, curvature, and gyrification measures)
- Quantitative mapping (eg. map your qT1 MRI data to a midthickness surface; extract laminar profiles perpendicular to this surface)

## Easy install

```bash
pixi global install hippunfold -c conda-forge -c khanlab -c bioconda

```
But see [Full Documentation](https://hippunfold.readthedocs.io/en/latest/?badge=latest) for details.

## What's new in version 2.0? 

- Better surface precision for capturing individual variability in gyral/sulcal/digitation patterning
- Retesselation of surfaces for more uniform face sizes
- Modified Laplace coordinate system to minimize distortion between folded and unfolded spaces
- Better robustness to topological breaks (a rare but annoying issue with previous version)
- New U-net models for neonates, Alzheimer's disease, contrast-agnostic processing, and more on the way
- Easier installation and instructions

An overview of these changes and epiricial testing coming soon to **biorxiv**

## Workflow 

![Pipeline Overview](https://raw.githubusercontent.com//khanlab/hippunfold/dev-v2.0.0/docs/images/HippUnfold_v2_overview.png)

For more information, see
**Full Documentation:**  [here](https://hippunfold.readthedocs.io/en/latest/?badge=latest)

## HippoMaps

Compare your data to a library of hippocampal surface maps from MRI, fMRI, histology, iEEG, and more using [HippoMaps](https://hippomaps.readthedocs.io/en/latest/). Also conisder contributing your own data to this open repo!

## Publications

### HippUnfold methods paper

- DeKraker, J., Haast, R. A., Yousif, M. D., Karat, B., Lau, J. C., Köhler, S., & Khan, A. R. (2022). Automated hippocampal unfolding for morphometry and subfield segmentation with HippUnfold. Elife, 11, e77945. [link](https://doi.org/10.7554/eLife.77945)
  -  **Please cite this if you use any version of HippUnfold)**
  
### Unfolded space registration and multihist7 atlas
- DeKraker Jordan, Palomero-Gallagher Nicola, Kedo Olga, Ladbon-Bernasconi Neda, Muenzing Sascha E.A., Axer Markus, Amunts Katrin, Khan Ali R., Bernhardt Boris, Evans Alan C. (2023) Evaluation of surface-based hippocampal registration using ground-truth subfield definitions eLife 12:RP88404 [link](https://doi.org/10.7554/eLife.88404.3)
  -  **Please cite this if you use HippUnfold version >= 1.3.0)**

### Commentary on surface-based hippocampal segmentation
- DeKraker J, Köhler S, Khan AR. Surface-based hippocampal subfield segmentation. Trends Neurosci. 2021 Nov;44(11):856-863. doi: 10.1016/j.tins.2021.06.005. Epub 2021 Jul 22. PMID: 34304910. [link](https://pubmed.ncbi.nlm.nih.gov/34304910/)

### Related papers

- DeKraker J, Ferko KM, Lau JC, Köhler S, Khan AR. Unfolding the hippocampus: An intrinsic coordinate system for subfield segmentations and quantitative mapping. Neuroimage. 2018 Feb 15;167:408-418. doi: 10.1016/j.neuroimage.2017.11.054. Epub 2017 Nov 23. PMID: 29175494. [link](https://pubmed.ncbi.nlm.nih.gov/29175494/)
- DeKraker J, Lau JC, Ferko KM, Khan AR, Köhler S. Hippocampal subfields revealed through unfolding and unsupervised clustering of laminar and morphological features in 3D BigBrain. Neuroimage. 2020 Feb 1;206:116328. doi: 10.1016/j.neuroimage.2019.116328. Epub 2019 Nov 1. PMID: 31682982. [link](https://pubmed.ncbi.nlm.nih.gov/31682982/)
- Karat BG, DeKraker J, Hussain U, Köhler S, Khan AR. Mapping the macrostructure and microstructure of the in vivo human hippocampus using diffusion MRI. Hum Brain Mapp. 2023 Nov;44(16):5485-5503. Epub 2023 Aug 24. PMID: 37615057; PMCID: PMC10543110.[link](https://doi.org/10.1002/hbm.26461)
- Ripart M, DeKraker J, Eriksson MH, Piper RJ, Gopinath S, Parasuram H, Mo J, Likeman M, Ciobotaru G, Sequeiros‐Peggs P, Hamandi K. Automated and interpretable detection of hippocampal sclerosis in temporal lobe epilepsy: AID‐HS. Annals of Neurology. 2025 Jan;97(1):62-75.
- Haast RA, Kashyap S, Ivanov D, Yousif MD, DeKraker J, Poser BA, Khan AR. Insights into hippocampal perfusion using high-resolution, multi-modal 7T MRI. Proceedings of the National Academy of Sciences. 2024 Mar 12;121(11):e2310044121.
- Eichert N, DeKraker J, Howard AF, Huszar IN, Zhu S, Sallet J, Miller KL, Mars RB, Jbabdi S, Bernhardt BC. Hippocampal connectivity patterns echo macroscale cortical evolution in the primate brain. Nature Communications. 2024 Jul 16;15(1):5963.
- Nichols ES, Karat BG, Grace M, Bezanson S, Khan AR, Duerden EG. Early life stress impairs hippocampal subfield myelination. Communications Biology. 2025 May 22;8(1):785.

### Critical tools 

- Huber LR, Poser BA, Bandettini PA, Arora K, Wagstyl K, Cho S, Goense J, Nothnagel N, Morgan AT, van den Hurk J, Müller AK. LayNii: A software suite for layer-fMRI. NeuroImage. 2021 Aug 15;237:118091.
- Isensee F, Jaeger PF, Kohl SA, Petersen J, Maier-Hein KH. nnU-Net: a self-configuring method for deep learning-based biomedical image segmentation. Nature methods. 2021 Feb;18(2):203-11.
- Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A. Sustainable data analysis with Snakemake. F1000Research. 2021 Apr 19;10:33.
