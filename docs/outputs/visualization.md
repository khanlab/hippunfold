# Visualization

## Freeview (volumes and surfaces)

[Freeview](https://surfer.nmr.mgh.harvard.edu/fswiki/FreeviewGuide/FreeviewGeneralUsage) is a powerful viewer that works well with both volume and surface data. It (currently) has several quirks worth noting, such as:
- All surfaces should be loaded before volumes (otherwise their orientation will be incorrect).
- Loading surface metric data can be difficult, see the examples below.
- Crashes can still occur occasionally, especially when trying to load surfaces or surface overlay data with the incorrect buttons.

### The example below will act as a guide to avoid common mistakes:
- open Freeview by simple typing `freeview` in the command line.
- Add a hippunfold surface:
  File > Load Surface > Navigate to a surface (eg. `hippunfold/sub-01/surf/sub-02_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii`).
- Add the correponding T1w image:
  File > Load Volume > Navigate to a volume (eg. `hippunfold/sub-01/anat/sub-02_desc-preproc_T1w.nii.gz`).
You should now see a hippocampal surface projected onto the coronal, sagittal, and axials views over a T1w image. You should also see a 3D model of the hippocampus. You can toggle visibility of each file by unticking or ticking it in the left panel. You can also adjust the ordering ov overlays here. 
- Add shape data as an overlay on the midthickness surface:
  With the surface file selected, use the left information panel to select Overlay > Load Generic and navigate to a metric file (eg. `hippunfold/sub-01/surf/sub-02_hemi-L_space-T1w_den-0p5mm_label-hipp_gyrification.shape.gii`).
You should now have a Configure button which you can use to adjust the windowing and colormap of this data.
- Subfield labels (`.label.gii` files) can be loaded similarly to metric data, but using the Annotation button on the left panel instead of Overlay.
- Add a segmentation image (eg. `hippunfold/sub-01/anat/sub-01_hemi-L_space-cropT1w_desc-subfield_dseg.nii.gz`). On the left panel, set Colormap > Lookup Table and then Select Lookup Table > Load Lookup Table > Navigate to `YOUR_HIPPUNFOLD_INSTALLATION_DIRECTORY/hippunfold/resources/desc-subfields_freeview_desg.tsv`. 
You should now see standardized subfield colours and names in the bottom panel when you mouse over a given subfield. 

### Other visualization tips and tricks:
- Consider adding dentate surfaces (`label-dentate`), unfolded surfaces (`space-unfolded`), or other overlay data (eg. `thickness`).
- Ensure surfaces and metric data are sampled with the same number of vertices (eg. `label-hipp` and den-'0p5mm'), and note that folded and unfolded surfaces will appear far apart in the 3D viewer and so you may need to zoom out quite far to navigate to them.
- Any metric data loaded on a folded (eg. `space-T1w`) surface can also be viewed on an unfolded surface that has the same number of vertices.
- Add Laplace coordinates (eg. `hippunfold/sub-01/coords/sub-01_dir-AP_hemi-L_space-cropT1w_label-hipp_desc-laplace_coords.nii.gz`) over your T1w image, then set the minimum windowing to 0.001 and tick the box `Clear background` in the left panel to maintain visibility of the T1w underlay.
- Note that `space-T1w` and `space-cropT1w` should appear in equivalent positions in Freeview, despite having a different field of view and resolution. If loading these data into Matlab or Python, you will need to first resample these images to the same reference image in order to index voxels from equivalent points. 

## HippUnfold Toolbox

The [HippUnfold Toolbox](https://github.com/jordandekraker/hippunfold_toolbox) provides examples and functions for mapping data onto hippocampal surfaces, plotting surfaces, and performing comparisons or statistical tests between subjects. Note that this can be done in other programs, like Connectome workbench, but these Python & Matlab tools should give an idea of how this can be done in a fully customizable fashion.


## ITK-SNAP (volumes)

[ITK-SNAP](http://www.itksnap.org/pmwiki/pmwiki.php) is a lightweight tool able to quickly open volumes, and is ideal for manual segmentation or edits. Segmentation images (`_dseg.nii.gz`) can be loaded as overlays and ITK-SNAP will create a 3D rendering of the contours of each label. Use the `space-cropT1w` or `space-cropT2w` images in the `anat/` folder to visualize one hemisphere at a time.

## Connectome Workbench (surfaces)

[Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) view tool allows for advanced visualization of surfaces and surface data. Loading one of the `.spec` files produced by HippUnfold into `wb_view` will allow you to visualize the hippocampus in native and unfolded configurations, and also overlay metric and label data (e.g. subfields and thickness).


