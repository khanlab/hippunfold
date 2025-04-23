# Output Files

The `PATH_TO_OUTPUT_DIR` folder contains a `logs` and `config` folder for troubleshooting, but for most purposes all the outputs of interest will be in subje t folders with the following structure:

    sub-1425/
    ├── anat
    ├── surf
    ├── metric
    ├── cifti
    └── qc


Briefly, `anat` contains preprocessed volumetric input images and output segmentations in nifti format, `surf` contains surface meshes in gifti format along with aggregate spec files, `metric` contains scalar metrics and labels in gifti format, `cifti` contains metrics aggregated over hemispheres, and `qc` contains snapshots and useful diagnostic information for quality control.

## anat

This folder contains input anatomical images that have been non-uniformity corrected,
motion-corrected, and, where appropriate, averaged and registered. 
In this example, a `T1w` image was used as a standard reference image, but a `T2w` was also registered and used in tissue segemntation:

    sub-001
     └── anat
         ├── sub-001_desc-preproc_T1w.nii.gz
         ├── sub-001_space-T1w_desc-preproc_T2w.nii.gz
         ├── sub-001_hemi-R_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz
         ├── sub-001_hemi-R_space-cropT1w_desc-preproc_T2w.nii.gz
         ├── sub-001_hemi-R_space-cropT1w_desc-subfields_atlas-multihist7_dseg.nii.gz
         └── sub-001_space-cropT1w_desc-subfields_atlas-multihist7_volumes.tsv

As per BIDS guidelines, `desc-preproc` refers to preprocessed input images, `space-T1w` refers to the volume to which the image is registered, `hemi` refers to the left or right hemisphere (only shown for the right in this example), and `dseg` (discrete-segmentation) images with `desc-subfields` contains subfield labels (coded as integers as described in the included `volumes.tsv` file). The subfield atlas used will also be included, by default as `atlas-bigbrain`. Note that HippUnfold does most intermediate processing in an unshown `space-corobl` which is cropped, upsampled, and rotated (available with `--output-space corobl`. Downsampling to the original `T1w` space can thus degrade the results and so they are also provided in a higher resolution `space-cropT1w` space which is ideal for conducting volumetry or morphometry measures with high precision and detail. 

For example, the following image shows a whole-brain `T1w` image, a
`space-cropT1w` overlay of the upsampled T2w image (centre square), and a similarly upsampled output
subfield segmentation (colour).

![image](../images/T1-T2-subfields_sag.png)


## surf

### surface meshes

Surface meshes (geometry files) are in `.surf.gii` format, and are
provided in both the anatomical space (`space=T1w, space=corobl`) and the unfolded space
(`space-unfold`). In each space, there are `inner`, `midthickness`,
and `outer` surfaces, which correspond to `white`, `midthickness`, and
`pial` for cortical surfaces:

    sub-{subject}
     └── surf
         └── sub-001_hemi-R_space-{T1w,unfold,corobl}_den-8k_{inner,midthickness,outer}.surf.gii

The following shows surfaces `inner`, `midthickness`, and `outer` in
yellow, orange, and red, respectively.

![image](../images/inner-mid-outer_sag.png)

### surface densities

Surfaces are provided in different density configurations, and are
labelled based on the approximate number of vertices in the hippocampus surface.  The default density is `8k`, which has an approximate vertex spacing of 0.5mm. Note that the dentate surface always has 1/4 the number of vertices, e.g. the dentate 8k surface actually has 2k vertices. 

There are also 
`2k` and `512` surfaces which have 1mm or 2mm spacing, respectively (suitable for
lower-resolution BOLD data). 

All surfaces of the same density (e.g. `8k`), in both
`space-T1w` and `space-unfold`, share the same mesh topology and have
corresponding vertices with each other. The vertex locations for
unfolded surfaces are identical for all subjects as well (note that this
of course is not the case for the anatomical `space=T1w,corobl` surfaces).

### label-dentate

HippUnfold v1.0.0 introduced `label-dentate` files which represent a distinct surface making up the dentate gyrus (reflecting its distinct topology from the rest of the cortex). 

These are illustrated in the following image (orange represents the usual hippocampal midthickness surface, while violet shows the new `dentate` surface):


### spec files

These files, along with the below metric and cifti files are packaged together for easy viewing in Connectome Workbench, `wb_view`, in a density-specific `.spec` file:

**TODO**: should remove space-T1w from spec file (isn't needed)..
    sub-{subject}
     └── surf
         └── sub-001_space-T1w_den-0p5mm_surfaces.spec  

## metrics

### shape metrics

In addition to the geometry files, surface-based shape metrics are
provided in `.shape.gii` format. The thickness, curvature and surface
area are computed using methods analogous to those used for cortical surfaces, based on
the surface geometry files, and are applicable to surfaces of the same density in any space.  The gyrification metric is the ratio of native to unfolded surface area, or
equivalently, the scaling or distortion factor when unfolding:

    sub-{subject}
     └── metric
         └── sub-001_hemi-{L,R}_den-8k_label-hipp_{thickness,curvature,gyrification}.shape.gii
         └── sub-001_hemi-{L,R}_den-8k_label-dentate_{curvature,gyrification}.shape.gii

These metrics are shown in both anatomical and unfolded space in the images
below. Note that these results are from group-averaged data and so
individual subject maps may show considerably more variability.

**TODO: update this image to show an individual example?

![image](../images/metrics.png)



### myelin maps

If your dataset has T1w and T2w images (and you are using `--modality=T1w` or `--modality=T2w`), then you can enable the generation of myelin maps as the ratio of T1w over T2w images. This division is done in the `corobl` space, and provides `myelin.shape.gii` surface metrics, and also includes these in the CIFTI and spec files. 

This option is enabled with the `--generate-myelin-maps` command-line option.

### coords

**TODO** add surf-based coords to the spec/target rules

Hippunfold also provides images that represent anatomical gradients
along the 3 principal axes of the hippocampus, longitudinal from
anterior to posterior, lamellar from proximal (dentate gyrus) to distal
(subiculum), and laminar from inner (SRLM) to outer. These are provided
in the metrics suffixed with `coords.shape.gii` with the direction indicated
by `dir-{direction}` as `AP`, `PD` or `IO`, and intensities from 0 to 1,
e.g. 0 representing the Anterior end and 1 the Posterior end.

**TODO**: update to show surface coords


### surface labels

The subfield labels from unfolded atlases are also provided for each 
subject, in `.label.gii` format. Analogous to the volume-based labels,
the name of the atlas (default: `multihist7`) is in the file name.

    sub-{subject}
     └── surf
         └── sub-001_hemi-{L,R}_den-8k_label-hipp_atlas-multihist7_subfields.label.gii


## cifti

In addition to lateralized `.shape.gii` and `.label.gii` metrics and labels,
we also provide data mapped from hemispheres in a single
file using the corresponding CIFTI formats, `.dscalar.nii` and `.dlabel.nii`. 

    sub-{subject}
     └── surf
         ├── sub-001_den-8k_label-{hipp,dentate}_{thickness,curvature,gyrification}.dscalar.nii
         └── sub-001_den-0p5mm_label-hipp_atlas-multihist7_subfields.dlabel.nii
**TODO** update this figure

![image](../images/dentate_cor.png)

## Additional Files

The top-level `PATH_TO_OUTPUT_DIR` contains additional folders:

    ├── logs
    ├── config
    └── .snakemake



The `config` folder contains the parameters as translated from the command-line interface, and passed along to the snakemake workflow by snakebids. 

The hidden `.snakemake` folder is created by snakemake, and contains a record of the code and parameters used, and paths to the inputs.

Workflow steps that write logs to file are stored in the `logs`
subfolder, with folders based on the rule name, and file names based on the wildcards (e.g. subject,
hemi, etc..).
