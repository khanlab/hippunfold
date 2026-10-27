# Quality Control

A complex pipeline can have multiple points of failure, so its imprtant to check the QC folder to ensure results make sense. This will contain (for example) the following files:

    hippunfold/
    └── sub-001
        └── qc
            ├── sub-001_desc-srcsink_qc.html
            ├── sub-001_desc-srcsink_qcchecks.tsv
            ├── sub-001_desc-surface_qc.html
            ├── sub-001_desc-surface_qcchecks.tsv
            ├── sub-001_from-T1w_to-CITI168_regqc.png
            ├── sub-001_hemi-L_label-hipp_desc-srcsink_qcmetrics.json
            ├── sub-001_hemi-L_label-hipp_desc-surface_qcmetrics.json
            ├── sub-001_hemi-L_desc-unetf3d_dice.tsv
            ├── sub-001_hemi-L_space-cropT1w_desc-subfields_atlas-histologyReference2023_dseg.png
            ├── sub-001_hemi-L_space-T1w_den-0p5mm_label-dentate_desc-subfields_midthickness.surf.png
            └── sub-001_hemi-L_space-T1w_den-0p5mm_label-hipp_desc-subfields_midthickness.surf.png

These files will be mirrored in the right (`hemi-R`) hemisphere, but we will skip that here for brevity. 

## Automated checks

A very fast quality check is to simply look at `sub-001_hemi-L_desc-unetf3d_dice.tsv` and ensure that this number is greater than `0.7`. Typically this criterion can be applied to all subjects, and subjects `<0.7` can simply be automatically discarded. A more in-depth explanation for this is provided at the bottom of this page. 

### QC reports and checks

Some processing steps record QC metrics and pass/fail checks as they run, in `*_qcmetrics.json` files. For each of these, an HTML report (`*_qc.html`) visualizes the metrics for the subject, and a summary of all checks is written to `*_qcchecks.tsv`. Open the HTML report in a web browser; it is a single self-contained file.

Checks have a severity of either *warning* or *error*. Failed warnings are highlighted in the report but do not stop processing. If an *error* check fails, HippUnfold still writes the report, then stops with an error naming the failed checks and the report to inspect. To finish processing anyway (e.g. after inspecting the report), re-run with `--qc_warn_only`.

Currently this includes:

- `desc-srcsink`: the anterior-posterior (AP) and proximal-distal (PD) source/sink terminals, which set the boundary conditions for the Laplace-Beltrami coordinates of each surface. They are fit to the outer boundary of the surface as four contiguous arcs, alternating AP and PD, so that each source is opposite its sink. The report shows the terminals along the boundary and on the surface, and the resulting AP and PD coordinates. It is an *error* if the terminals are in the wrong order for the hemisphere (swapped terminals or a flipped surface), and a *warning* if the terminals from the segmentation are fragmented, had to be changed substantially, or are far from the surface boundary.
- `desc-surface`: the topology and mesh quality of each native surface. The hippocampal (and dentate) surface should be topologically a disk: one connected, orientable piece with genus *g* = 0 (no handles) and *b* = 1 boundary loop, so its Euler characteristic is χ = V − E + F = 2 − 2*g* − *b* = 1. The report gives χ, *g* and *b* for each surface, and lists every defect found, numbered and located on the surface (position in mm, and AP/PD coordinates): holes, non-manifold or misoriented elements, self-intersections of the midthickness, inner or outer surfaces, degenerate or sliver triangles, interior edges joining two boundary vertices (which make faces collapse when the terminal arcs of the boundary are fixed), and faces that are folded over or collapsed to zero area in the unfolded (AP, PD) mapping. The face and vertex indices of each defect are in the `*_qcmetrics.json` files. It is an *error* if a surface has more than one connected component, no boundary, or is non-orientable; all other defects are *warnings*, reported whenever there are any (the thresholds can be changed in `surface_qc_thresholds` in `rules/qc.smk`).

  The midthickness is also checked at intermediate stages of surface processing (`stage-isosurf`: the raw isosurface, before holes are filled and defects repaired; `stage-repaired`: after repair and decimation, before the triangle quality is improved by edge flips and tangential relaxation), in `*_desc-surface_stage-*_qcmetrics.json`. The report compares each stage with the final surface, with a defect map for each. Only the final surface is used for the pass/fail checks.

To compare the QC of two runs, e.g. before and after a change to the workflow, use:

    python hippunfold/workflow/scripts/qc_compare_runs.py <before_dir> <after_dir> --qc surface --labels before after -o comparison.html

Defect maps are included when the native surfaces were kept (with `--notemp`).

To better understand some common points of failure, let's walk through the other files contained here:

## Catch-all failures

A good way to be sure that nothing has failed is also to look at the final results, as in the following snapshot images.

`sub-001_hemi-L_space-cropT1w_desc-subfields_atlas-histologyReference2023_dseg.png`:

![image](../images/QC/sub-001_hemi-L_space-cropT1w_desc-subfields_atlas-histologyReference2023_dseg.png)

Here we see a volumetric output of subfield parcellations. If you have some hippocampal expertise then you may be able to see right away whether the results make sense. Some common errors include missed parts of the hippocampus and, sometimes, parts of collateral sulcus are included as part of the hippocapus. This is because of gross or "catastrophic" failures in the UNet tissue classification, which should be caught be the above automated check. Still, this collage of sample slices should give a good sense for whether/how the segmentation has failed. 

`sub-001_hemi-L_space-T1w_den-0p5mm_label-hipp_desc-subfields_midthickness.surf.png`:

![image](../images/QC/sub-001_hemi-L_space-T1w_den-0p5mm_label-hipp_desc-subfields_midthickness.surf.png)

`sub-001_hemi-L_space-T1w_den-0p5mm_label-dentate_desc-subfields_midthickness.surf.png`:

![image](../images/QC/sub-001_hemi-L_space-T1w_den-0p5mm_label-dentate_desc-subfields_midthickness.surf.png)

In these two images we see the output hippocampal and dentate gyrus surfaces, which should roughly resemble those shown in the HippUnfold manuscript (though they may have some differences in shape and/or sulcal/gyral patterning due to inter-individual variance!). 

Sometimes the dentate gyrus surface may show a few large triangles that obscure the rest of the surface, as in the example here. This typically arises when just a few vertices are badly misplaced, and so typically this doesn't have much impact on quantitative results but it can look quite ugly in a figure. This can arise because the structure is small and contains few voxels, none of which have a laplace field close to a given value and so the xyz coordinate of a given laplace value is extrapolated badly. We will try to fix this issue in future releases.

## Registration/cropping failures

Though it seems simple, just getting the image roughly aligned to a template so we can then crop around the left and right hippocampi for subsequent processing (in `space-corobl`) can sometimes fail. We can examine this in the file `sub-001_from-T1w_to-CITI168_regqc.png`:

![image](../images/QC/sub-001_from-T1w_to-CITI168_regqc.png)

In this image, the subject tissue boundaries (e.g. grey matter, white matter, CSF, skull) are overlaid on a standard template. This is a linear alignment, so while gyri/sulci may not be well aligned, we should still see fairly good overlap that is sufficient to derive the position of the hippocampus in the image. This most often fails if the input images are in an unexpected orientation (e.g. this sometimes arises in ex-vivo scanning), and can sometimes be imporved by first running FSL's [reorient2std](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils). 

## Checking the automated check

The `Dice > 0.7` rule is a good automatic check because HippUnfold typically shows high sensitivity to details in hippocampal structure, but when it fails it often does so catastrophically. These failures are very obvious, and can often be shown by this simple Dice score that compares a whole mask of the hippocampus to a whole mask of a hippocampus generated by a more conventional but less precise measure: diffeomorphic registration to a standard template. The Dice score should be good between the two methods, but will never be close to perfect (`>0.9`) since the HippUnfold method will differ as its more precise. 

This issue typically only occurs in ~1% of cases. 
