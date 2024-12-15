surf_thresholds={'inner': 0.05, 'outer':0.95, 'midthickness':0.5}


rule all_nativesurf:
    input:
        thickness_gii=expand(bids(
            root=root,
            datatype="surf_",
            suffix="thickness.shape.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),subject='1425',desc='equivol',hemi='R',label='hipp')
      


rule get_hipp_mask_for_meshing:
    input:
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=' '.join([str(lbl) for lbl in config["laplace_labels"]["AP"]["gm"]])

    output:
        mask=bids(
            root=root,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    shell:
        'c3d {input} -retain-labels {params} -binarize {output}' 



rule gen_native_hipp_mesh:
    input:
        nii=bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="hipp",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        mask=bids(
            root=root,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    params:
        threshold=lambda wildcards: surf_thresholds[wildcards.surfname],
        decimate_percent=0 # not currently working right now
    output:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname,midthickness}.surf.gii",
            space="corobl",
            desc="{desc}nostruct",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    #container (will need pyvista dependency)
    script:
        "../scripts/gen_isosurface.py"

rule update_native_hipp_mesh_structure:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}nostruct",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname,midthickness}.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "cp {input} {output} && wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule warp_native_mesh_to_unfold:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        warp_native2unfold=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="{label}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="surface"
        ),
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="FLAT",
    output:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="unfold",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp_native2unfold} {output.surf_gii} && "
        "wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule resample_subfields_to_native_surf:
    input:
        label_gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="subfields.label.gii",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_hipp",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="midthickness.surf.gii",
            space="unfold",
            desc="{desc}",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf_",
            fromdensity="{density}",
            desc="{desc}",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container: None
#        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -label-resample {input.label_gii} {input.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check'

rule resample_native_surf:
    input:
        native_corobl=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surf_name}.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_hipp",
            "tpl-avg_space-unfold_den-{density}_{surf_name}.surf.gii",
        ),
        native_unfold=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surf_name}.surf.gii",
            space="unfold",
            desc="{desc}",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    output:
        native_corobl_resampled=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surf_name}.surf.gii",
            space="corobl",
            den="{density}",
            desc="{desc}",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    container: None
#        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -surface-resample {input.native_corobl} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_corobl_resampled} -bypass-sphere-check'



rule compute_halfthick_mask:
    input:
        coords=bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="hipp",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        mask=bids(
            root=root,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    params:
        threshold_tofrom=lambda wildcards: "0.5 1" if wildcards.inout == 'inner' else "0 0.5"
    output:
        nii=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="mask.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    shell:
        "c3d {input.coords} -threshold {params.threshold_tofrom} 1 0 {input.mask} -multiply -o {output}"

rule register_midthickness:
    input:
        fixed=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="mask.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        moving=bids(
            root=root,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    output:
        warp=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="xfm.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),       
    shell:
        'greedy -d 3 -i {input.fixed} {input.moving} -n 30 -o {output.warp}'

    

rule apply_halfsurf_warp_to_img:
    input:
        fixed=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="mask.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        moving=bids(
            root=root,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
        warp=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="xfm.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),       
    output:
        warped=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="warpedmask.nii.gz",
            desc="{desc}",
            to_="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    shell:
        'greedy -d 3  -rf {input.fixed} -rm {input.moving} {output.warped}  -r {input.warp} '


rule convert_warp_from_itk_to_world:
    input:
        warp=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="xfm.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),       

    output:
        warp=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="xfmras.nii.gz",
            desc="{desc}",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),       
    shell:
        'wb_command -convert-warpfield -from-itk {input} -to-world {output}'


rule warp_midthickness_to_inout:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="midthickness.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        warp=bids(
            root=root,
            datatype="surf_",
            dir="IO",
            label="hipp",
            suffix="xfmras.nii.gz",
            desc="{desc}",
            to="{surfname}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),       
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname,inner|outer}.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp} {output.surf_gii} && "
        "wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"




rule calc_thickness:
    input:
        inner=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="inner.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        outer=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="outer.surf.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf_",
            suffix="thickness.shape.gii",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"

 
