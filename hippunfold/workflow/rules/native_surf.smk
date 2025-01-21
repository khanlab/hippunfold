# This is one big rules file for now, but could split it up according to the
# dividers ( --- ) that describe each conceptual part of the pipeline

# --- parameters - will put these in config eventually

surf_thresholds = {"inner": 0, "outer": 1, "midthickness": 0.5}


# this is for the mapping from inner to outer
gm_labels = {
    "hipp": config["laplace_labels"]["IO"]["gm"],
    "dentate": config["laplace_labels"]["PD"]["sink"],
}

# appends the coords with these regions set to +1.1 for the meshing
sink_labels = {"hipp": config["laplace_labels"]["IO"]["sink"], "dentate": [2]}

# appends the coords with these regions set to +1.1 for the meshing
src_labels = {"hipp": config["laplace_labels"]["IO"]["src"], "dentate": [1]}

# sets these to nan in the coords for the meshing
nan_labels = {
    "hipp": config["laplace_labels"]["AP"]["sink"]
    + config["laplace_labels"]["AP"]["src"]
    + config["laplace_labels"]["PD"]["sink"]
    + config["laplace_labels"]["PD"]["src"],
    "dentate": [
        0
    ],  # TODO: this requires labels we don't produce yet -- namely, those at the boundary between  PDcoord~0.9-1  and SRLM, and between PDcoord~0.9-1 and BG
}

desc_io = {
    "hipp": "equivol" if "equivolume" in config["laminar_coords_method"] else "laplace",
    "dentate": "laplace",
}

unfoldreg_method = "greedy"  # choices: ["greedy","SyN"]

unfoldreg_padding = "64x64x0vox"


ruleorder: resample_native_surf_to_std_density > cp_template_to_unfold
ruleorder: atlas_label_to_unfold_nii > atlas_metric_to_unfold_nii


# --- isosurface generation ---


rule get_label_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in gm_labels[wildcards.label]]
        ),
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                suffix="mask.nii.gz",
                space="corobl",
                desc="GM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -retain-labels {params} -binarize {output}"


rule get_sink_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in sink_labels[wildcards.label]]
        ),
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                suffix="mask.nii.gz",
                space="corobl",
                desc="sink",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_src_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in src_labels[wildcards.label]]
        ),
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                suffix="mask.nii.gz",
                space="corobl",
                desc="src",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_nan_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in nan_labels[wildcards.label]]
        ),
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                suffix="mask.nii.gz",
                space="corobl",
                desc="nan",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule gen_native_mesh:
    input:
        coords=lambda wildcards: bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="coords.nii.gz",
            desc=desc_io[wildcards.label],
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        nan_mask=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="nan",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        sink_mask=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="sink",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        src_mask=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="src",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        threshold=lambda wildcards: surf_thresholds[wildcards.surfname],
        decimate_percent=0,  # not enabled
    output:
        surf_gii=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="{surfname,midthickness}.surf.gii",
                space="corobl",
                desc="nostruct",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/gen_isosurface.py"


rule update_native_mesh_structure:
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="nostruct",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname,midthickness}.surf.gii",
            space="corobl",
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


rule smooth_surface:
    """ slight smoothing of surface to improve curvature estimation """
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        smoothing_strength=0.8,
        smoothing_iterations=10,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="smoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-smoothing {input} {params.smoothing_strength} {params.smoothing_iterations} {output}"


# --- creating unfold surface from native anatomical, including post-processing


rule laplace_beltrami:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        seg=get_labels_for_laplace,
    params:
        src_labels=lambda wildcards: config["laplace_labels"],
    output:
        coords_AP=bids(
            root=work,
            datatype="coords",
            dir="AP",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        coords_PD=bids(
            root=work,
            datatype="coords",
            dir="PD",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/laplace-beltrami.py"


def get_unfold_z_level(wildcards):
    extent = float(config["unfold_vol_ref"][wildcards.label]["extent"][-1])
    return surf_thresholds[wildcards.surfname] * extent


rule warp_native_mesh_to_unfold:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        coords_AP=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="{label}",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        coords_PD=bids(
            root=work,
            datatype="coords",
            dir="PD",
            label="{label}",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    params:
        z_level=get_unfold_z_level,
        vertspace=lambda wildcards: config["unfold_vol_ref"][wildcards.label],
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldraw",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    script:
        "../scripts/rewrite_vertices_to_flat.py"


rule update_unfold_mesh_structure:
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldraw",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="FLAT",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
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


rule heavy_smooth_unfold_surf:
    """ this irons out the surface to result in more even
        vertex spacing. the resulting shape will be more 
        individual (e.g. the surface area in unfolded space 
        would be similar to native) -- TODO: maybe this is a good 
        way to determine smoothing strenghth and iterations, e.g. 
        use the surface area and vertex spacing as objectives.."""
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        strength=0.1,
        iterations=1000,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldsmoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-smoothing {input} {params.strength} {params.iterations} {output}"


# --- creating inner/outer surfaces from native anatomical (using 3d label deformable registration)


rule compute_halfthick_mask:
    input:
        coords=lambda wildcards: bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="coords.nii.gz",
            desc=desc_io[wildcards.label],
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        mask=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        threshold_tofrom=lambda wildcards: "0.5 1"
        if wildcards.inout == "inner"
        else "0 0.5",
    output:
        nii=temp(
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="mask.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.coords} -threshold {params.threshold_tofrom} 1 0 {input.mask} -multiply -o {output}"


rule register_midthickness:
    input:
        fixed=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        moving=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        warp=temp(
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="xfm.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "greedy -d 3 -i {input.fixed} {input.moving} -n 30 -o {output.warp}"


rule apply_halfsurf_warp_to_img:
    input:
        fixed=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        moving=bids(
            root=work,
            datatype="anat",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        warp=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="xfm.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    output:
        warped=temp(
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="warpedmask.nii.gz",
                to_="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "greedy -d 3  -rf {input.fixed} -rm {input.moving} {output.warped}  -r {input.warp} "


# TODO: rename the warps here according to existing custom (type=itk -> type=ras)
rule convert_inout_warp_from_itk_to_world:
    input:
        warp=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="xfm.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    output:
        warp=temp(
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="xfmras.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule warp_midthickness_to_inout:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        warp=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="xfmras.nii.gz",
            to="{surfname}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname,inner|outer}.surf.gii",
            space="corobl",
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


# --- affine transforming anatomical surfaces from corobl to other (T1w, T2w) spaces


# warp native surface from corobl to T1w/T2w
rule affine_gii_corobl_to_modality:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="ras"
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="{native_modality}",
            hemi="{hemi}",
            label="{autotop,hipp|dentate}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}"


# --- calculating metrics on native anatomical surfaces


rule calculate_surface_area:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule calculate_legacy_gyrification:
    """new gyrification is ratio of nativearea to unfoldarea (e.g. surface scaling or distortion factor.
    this should be proportional by a constant, to the earlier gyrification on 32k surfaces."""
    input:
        native_surfarea=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        unfold_surfarea=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="gyrification.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea}"


rule calculate_curvature:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            desc="smoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="curvature.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-curvature {input} -mean {output}"


rule calculate_thickness:
    input:
        inner=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        outer=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="thickness.shape.gii",
            space="corobl",
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


# --- unfolded registration


# for unfold_reg, need to get data into 2d image - Jordan's previous code used metric-to-volume-mapping with unfoldiso inner and outer surfaces
# we can use same approach, just calculating the metrics on native, and using the unfold inner/outer as ribbon


rule pad_unfold_ref:
    """pads the unfolded ref in XY (to improve registration by ensuring zero-displacement at boundaries)
        and to deal with input vertices that are just slightly outside the original bounding box. 
       The ribbon-constrained method will be used to resample surface metrics into this reference, then
        to get a 2D slice we will take the central slice (the next rule creates the single-slice reference
        for this)."""
    input:
        ref_nii=bids(
            root=work,
            space="unfold",
            label="{label}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
    params:
        padding=f"-pad {unfoldreg_padding} {unfoldreg_padding}",
    output:
        ref_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="padded",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.ref_nii} -scale 0 -shift 1 {params.padding} 0 -o {output.ref_nii}"


rule extract_unfold_ref_slice:
    """This gets the central-most slice of the unfold volume, for obtaining a 2D slice"""
    input:
        ref_3d_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="padded",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        ref_2d_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.ref_3d_nii} -slice z 50% -o {output.ref_2d_nii}"


rule native_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs"""
    input:
        metric_gii=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        ref_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        interp="-nearest-vertex 0.3",
    output:
        metric_nii=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.midthickness_surf} {input.ref_nii} {output.metric_nii} "
        " -ribbon-constrained {input.inner_surf} {input.outer_surf}"


rule atlas_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs. 
        This rule is for the surface template"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        metric_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["metric_gii"].format(**wildcards),
        inner_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="inner", **wildcards
        ),
        outer_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="outer", **wildcards
        ),
        midthickness_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        metric_nii=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -metric-to-volume-mapping {params.metric_gii} {params.midthickness_surf} {input.ref_nii} {output.metric_nii} "
        " -ribbon-constrained {params.inner_surf} {params.outer_surf}"


def get_fixed_images_unfoldreg(wildcards):
    unfoldreg_metrics = config["atlas_files"][wildcards.atlas]["metric_wildcards"]

    return expand(
        bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        metric=unfoldreg_metrics,
        **wildcards
    )


def get_moving_images_unfoldreg(wildcards):
    unfoldreg_metrics = config["atlas_files"][wildcards.atlas]["metric_wildcards"]
    return expand(
        bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
        metric=unfoldreg_metrics,
        **wildcards
    )


rule unfoldreg_antsquick:
    """performs 2D registration in unfolded space as in [reference paper]
    this is done in a shadow directory to get rid of the tmp files generated by ants.

    Note: fixed and moving are swapped as compared to v1 unfoldreg.

    TODO: this currently uses NMI as a metric, which doesn't work very well for
        deformable registraiton (won't warp very much).. should also use the 
        antsRegistration tool directly instead of the simple wrapper to provide
        more flexibility with specifying parameters. """
    input:
        fixed_images=get_fixed_images_unfoldreg,
        moving_images=get_moving_images_unfoldreg,
    params:
        antsparams="-d 2 -t so -o tmp",
        fixed_args=lambda wildcards, input: " ".join(
            ["-f {img}".format(img=img) for img in input.fixed_images]
        ),
        moving_args=lambda wildcards, input: " ".join(
            ["-m {img}".format(img=img) for img in input.moving_images]
        ),
        cmd_copy_warps=lambda wildcards, output: " && ".join(
            [
                f"cp tmp1Warp.nii.gz {output.warp}",
                f"cp tmp1InverseWarp.nii.gz {output.invwarp}",
            ]
        ),
    output:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="{atlas}",
            to="native",
            space="unfold",
            type_="itk2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="native",
            to="{atlas}",
            space="unfold",
            type_="itk2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="unfoldreg.txt",
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    shadow:
        "minimal"
    threads: 16
    resources:
        mem_mb=16000,
        time=10,
    shell:
        "antsRegistrationSyNQuick.sh {params.antsparams} {params.fixed_args} {params.moving_args} &> {log} && "
        "{params.cmd_copy_warps}"


rule unfoldreg_greedy:
    """performs 2D registration in unfolded space using greedy.
    Greedy produces the forward map (can optionally produce inverse maps)
    so is not symmetric, but can usually produce more precise warps
    if symmetry isn't desired"""
    input:
        fixed_images=get_fixed_images_unfoldreg,
        moving_images=get_moving_images_unfoldreg,
    params:
        in_images=lambda wildcards, input: " ".join(
            [
                f"-i {fix} {mov}"
                for fix, mov in zip(input.fixed_images, input.moving_images)
            ]
        ),
        metric="-m NCC 4x4x4",  #4x4x4 implies a 9x9x9 window (ie 4 on either side)
        regularization="-s 1.732vox 0.707vox",  #gradient sigma, warp sigma (defaults: 1.732vox, 0.707vox)
    output:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="{atlas}",
            to="native",
            space="unfold",
            type_="itk2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="unfoldreg_greedy.txt",
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    shadow:
        "minimal"
    threads: 16
    resources:
        mem_mb=16000,
        time=10,
    shell:
        "greedy -d 2 {params.metric} {params.regularization} {params.in_images} -o {output.warp}"


rule extend_warp_2d_to_3d:
    """ we have a 2d warp, in space of innersurf 2d slice, 254x126. we want to 
        1) make it a 3d warp (create a z-displacement, as zeros)
        2) extend the 2d warp across z slices, (so each z slice is transformed identically)
        3) fill in the central part of the array, when going from 254x126 to 256x128

        Note that this rule is hardcoded for 256x128xN ref and 254x126 2d warp; 
    """
    input:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="{desc}",
            from_="{from}",
            to="{to}",
            space="{space}",
            type_="itk2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="padded",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="{desc}",
            from_="{from}",
            to="{to}",
            space="{space}",
            type_="itk",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    script:
        "../scripts/convert_warp_2d_to_3d.py"


rule convert_unfoldreg_warp_from_itk_to_world:
    input:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="{desc}",
            from_="{to}",
            to="{from}",
            space="{space}",
            type_="itk",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="{desc}",
            from_="{from}",
            to="{to}",
            space="{space}",
            type_="surface",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


def get_unfold_ref_name(wildcards):
    if (
        wildcards.label in config["atlas_files"][config["atlas"]]["label_wildcards"]
        and config["no_unfolded_reg"] == False
    ):
        return "unfoldreg"
    else:
        return "unfold"


def get_unfold_ref(wildcards):
    """function to return either unfoldreg or unfold ref mesh, depending on whether
    unfoldreg can be performed (based on atlas wildcards)"""

    return bids(
        root=root,
        datatype="surf",
        suffix="midthickness.surf.gii",
        space=get_unfold_ref_name(wildcards),
        hemi="{hemi}",
        label="{label}",
        **inputs.subj_wildcards
    )


rule warp_unfold_native_to_unfoldreg:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc=unfoldreg_method,
            from_="native",
            to=config["atlas"],
            space="unfold",
            type_="surface",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="FLAT",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldreg",
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


# --- resampling using the unfoldreg surface to (legacy) standard densities (0p5mm, 1mm, 2mm, unfoldiso)
rule resample_atlas_subfields_to_std_density:
    """ resamples subfields from a custom atlas mesh (e.g. from an atlas anatomical/native mesh
    warped to unfold space) to a standard density"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        atlas_unfold=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -label-resample {params.label_gii} {params.atlas_unfold} {input.ref_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


rule resample_native_surf_to_std_density:
    input:
        native=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        native_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="{space,unfoldreg|corobl}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-resample {input.native} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check"


rule resample_native_metric_to_std_density:
    input:
        native_metric=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.{metrictype}.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        metric_resampled=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.{metrictype,shape|func}.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -metric-resample {input.native_metric} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.metric_resampled} -bypass-sphere-check"


rule cp_surf_to_root:
    input:
        native_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        native_resampled=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


# --- resampling from atlasnative to native vertices
rule resample_atlas_subfields_to_native_surf:
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        native_unfold=get_unfold_ref,
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        ref_unfold=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label,hipp}",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -label-resample {params.label_gii} {params.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


# --- rule for converting atlasnative subfields to unfold nii (e.g. analogous to unfoldiso standard space)
rule atlas_label_to_unfold_nii:
    """converts metric .gii files to .nii - this is a band-aid fix since we make use of the volumetric
        subfield labels in unfold space for painting the native volumetric dseg. Ie this uses the unfold volumetric (or effectively unfoldiso)
        to perform mapping from atlas to subject. better approach would be to adjust the downstream  volumetric dseg function to 
        make use of gifti labels instead.. 
"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_nii=bids(
            root=work,
            space="unfold",
            label="{label}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        inner_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="inner", **wildcards
        ),
        outer_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="outer", **wildcards
        ),
        midthickness_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {params.label_gii} {params.midthickness_surf} {input.ref_nii} {output.label_nii} "
        " -nearest-vertex 1000"


# --- cifti for native surfs (could use generic rule if we include den wildcard)


def get_inputs_cifti_metric_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_metric"] = (
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="L",
                label="{autotop}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_metric"] = (
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="R",
                label="{autotop}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dscalar_metric_cifti_native:
    input:
        unpack(get_inputs_cifti_metric_native),
    params:
        cmd=get_cmd_cifti_metric,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.dscalar.nii",
            space="{space}",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


def get_inputs_cifti_label_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_label"] = (
            bids(
                root=root,
                datatype="surf",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="L",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_label"] = (
            bids(
                root=root,
                datatype="surf",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="R",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dlabel_cifti_subfields_native:
    input:
        unpack(get_inputs_cifti_label_native),
    params:
        cmd=get_cmd_cifti_label,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
            space="{space}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


# --- create spec file (adapted from rule in gifti.smk)


rule create_spec_file_hipp_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        subfields=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **inputs.subj_wildcards
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        cifti_labels=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="subfields.dlabel.nii",
                atlas="{atlas}",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,hipp}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmds}"


rule create_spec_file_dentate_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,dentate}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmds}"


rule merge_lr_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                suffix="surfaces.spec",
                hemi="{hemi}",
                space="{space}",
                label="{autotop}",
                **inputs.subj_wildcards
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            space="{space}",
            suffix="surfaces.spec",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


rule merge_hipp_dentate_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                suffix="surfaces.spec",
                space="{space}",
                label="{autotop}",
                **inputs.subj_wildcards
            ),
            autotop=config["autotop_labels"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            space="{space}",
            suffix="surfaces.spec",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


rule cp_native_surf_to_root:
    input:
        native=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        native=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"
