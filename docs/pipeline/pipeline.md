# Pipeline Details


This section describes the HippUnfold workflow, that is, 
the steps taken to produce the intermediate and final files.
HippUnfold is a Snakemake workflow, and thus the workflow is a 
 directed acyclic graph (DAG) that is automatically configured based
on a set of rules. 

Below is a example *simplified* visualization of the workflow DAG
for the `--modality T1w` workflow.  Each rounded rectangle in the DAG 
represents a *rule*, that is, some code or script that produces 
output file(s), and the arrows represent file inputs and outputs to these 
rules.  It is *simplified* in that multiple 
instances of each rule are not shown, e.g. the `run_inference` rule 
runs on both left and right hemispheres (`hemi=L`, `hemi=R`),
but only one `run_inference` box is shown here.  

<img src="../../hippunfold/dags/out_rulegraph/T1w.svg" width="800px">

Although it may still look very complex (click on the image to 
enlarge), it is also  organized into groups of rules, 
each representing the main phases of the workflow. Each
grouped set of rules also exist in separate rule files,
 which can be found in
 the [rules sub-folder](http://github.com/khanlab/hippunfold/tree/master/hippunfold/workflow/rules)  
in the workflow source. For example, the [preproc_t1](http://github.com/khanlab/hippunfold/tree/master/hippunfold/workflow/rules/preproc_t1.smk)  file contains the rules related to pre-processing the T1w images, and these are 
grouped together in the above diagram by a blue rectangle labelled `preproc_t1`. 


The main phases of the workflow are described in the sections below, zooming in
on the rules used in each blue rectangle, one at a time.


## Pre-processing

The pre-processing workflow for HippUnfold is generated based on the input data (e.g. whether 
there are multiple T2w images or a single T2w image), what modality is used
 (e.g. `--modality T1w` or `--modality T2w`), and what optional arguments are
 specified (e.g. `--t1-reg-template`). 

## T1w pre-processing

T1w images are imported, intensity-corrected using N4, and linearly registered 
to the template image (default: CITI168 - an HCP T1w template). An existing transformation to 
align the images in a coronal oblique (`space-corobl`) orientation is concatenated, and 
this space is used to define the left and right hippocampus bounding boxes in 0.3mm isotropic space. The left 
hippocampus subvolume is left-right flipped at this stage too (subsequent steps in the `corobl` space operate
on both the `hemi-R` and `hemi-Lflip` images).

![T1w workflow preproc t1](../../hippunfold/dags/out_dag/T1w.preproc_t1.svg)


## T2w pre-processing

T2w images are processed similarly, except the T2w version of the template is used. If multiple T2w images
exist, these are motion-corrected and averaged prior to N4 correction. The diagram below shows the T2w pre-
processing workflow for a dataset with three T2w runs.

![T2w multi preproc_t2](../../hippunfold/dags/out_dag/T2w_multi.preproc_t2.svg)


For T2w images where template registration is failing (e.g. because the T2w images have a limited FOV),
the `--t1-reg-template` option can be used, and will perform template registration with the T1w images, along with 
a within-subject registration of the T2w to the T1w, concatenating all the transforms. This is shown in the diagrams below (with a single T2w image in this case):

![T2w t1-reg-template preproc_t1](../../hippunfold/dags/out_dag/T2w_t1-reg-template.preproc_t1.svg)
![T2w t1-reg-template preproc_t2](../../hippunfold/dags/out_dag/T2w_t1-reg-template.preproc_t2.svg)


## U-net segmentation

![U-net workflow](../../hippunfold/dags/out_dag/T1w.nnunet.svg)

![U-net architecture](../images/nnUnet_hippunfold.png)

## Template-based shape injection


The following diagram shows the workflow, but simplified to contain one hemisphere (`--hemi R`), and excluding the dentate gyrus.

![Shape injection](../../hippunfold/dags/out_dag/T1w_hemi-R_hipponly.shape_inject.svg)


## Laplace & equivolume coordinates

The following diagram shows the workflow, but simplified to contain one hemisphere (`--hemi R`), and excluding the dentate gyrus.

![Laplace coordinates](../../hippunfold/dags/out_dag/T1w_hemi-R_hipponly.autotop.svg)

## Subfields processing

The following diagram shows the workflow, but simplified to contain one hemisphere (`--hemi R`), and excluding the dentate gyrus.

![Subfields](../../hippunfold/dags/out_dag/T1w_hemi-R_hipponly.subfields.svg)


## Generating warp files

The following diagram shows the workflow, but simplified to contain one hemisphere (`--hemi R`), and excluding the dentate gyrus.

![Generating warps](../../hippunfold/dags/out_dag/T1w_hemi-R_hipponly.warps.svg)

## Surface processing

The following diagram shows the workflow, but simplified to contain one hemisphere (`--hemi R`), and excluding the dentate gyrus.

![Surface processing](../../hippunfold/dags/out_dag/T1w_hemi-R_hipponly.gifti.svg)


## Additional steps

Resampling to output resolution, quality control snapshot generation, and archiving the work folder are steps 
that are also carried out by the workflow, but the DAGs are now shown here because of the many inputs/outputs, and the linear
workflow structure.






