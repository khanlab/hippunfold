# This is one big rules file for now, but could split it up according to the
# dividers ( --- ) that describe each conceptual part of the pipeline

# --- parameters - will put these in config eventually


wildcard_constraints:
    label="[0-9a-zA-Z]+",
    metric="[0-9a-zA-Z]+",
    native_modality="[0-9a-zA-Z]+",


surf_thresholds = {"inner": 0, "outer": 1, "midthickness": 0.5}

# this is for the mapping from inner to outer
gm_labels = {
    "hipp": config["laplace_labels"]["IO"]["gm"],
    "dentate": config["laplace_labels"]["PD"]["sink"],
}

desc_io = {
    "hipp": "equivol" if "equivolume" in config["laminar_coords_method"] else "laplace",
    "dentate": "laplace",
}

resample_from_density = (
    "unfoldiso"  # what density to resample metrics (ie subfields) from
)


# could allow user to select additional metrics (e.g. from bold, dwi ...)
unfoldreg_metrics = ["thickness", "curvature", "gyrification"]


ruleorder: resample_native_surf_to_std_density > cp_template_to_unfold
ruleorder: unfold_spring_model > warp_unfold_native_to_unfoldreg
ruleorder: atlas_label_to_unfold_nii > atlas_metric_to_unfold_nii  #temporary until we change from subfields.nii to desc-subfields_dseg.nii for unfold space   


# --- isosurface generation ---


rule prep_hipp_coords_for_meshing:
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
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: config["laplace_labels"]["IO"]["gm"],
        gm_noDG_labels=lambda wildcards: config["laplace_labels"]["IO"]["gm_noDG"],
        src_labels=lambda wildcards: config["laplace_labels"]["IO"]["src"],
        sink_labels=lambda wildcards: config["laplace_labels"]["IO"]["sink"],
    output:
        coords=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="coords.nii.gz",
                space="corobl",
                desc="formesh",
                hemi="{hemi}",
                label="{label,hipp}",
                **inputs.subj_wildcards
            )
        ),
        mask=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="mask.nii.gz",
                space="corobl",
                desc="formesh",
                hemi="{hemi}",
                label="{label,hipp}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/prep_hipp_coords_for_meshing.py"


rule prep_dentate_coords_for_meshing:
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
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: config["laplace_labels"]["IO"]["gm"],
    output:
        coords=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="coords.nii.gz",
                space="corobl",
                desc="formesh",
                hemi="{hemi}",
                label="{label,dentate}",
                **inputs.subj_wildcards
            )
        ),
        mask=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="mask.nii.gz",
                space="corobl",
                desc="formesh",
                hemi="{hemi}",
                label="{label,dentate}",
                **inputs.subj_wildcards
            )
        ),
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/prep_dentate_coords_for_meshing.py"


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
    shell:
        "c3d {input} -retain-labels {params} -binarize {output}"


rule gen_native_mesh:
    input:
        nii=lambda wildcards: bids(
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
    # container (will need pyvista dependency)
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
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
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


# --- creating unfold surface from native anatomical


rule warp_native_mesh_to_unfold:
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
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
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp_native2unfold} {output.surf_gii} && "
        "wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


# --- creating inner/outer surfaces from native anatomical (using 3d label deformable registration)


rule compute_halfthick_mask:
    input:
        coords=bids(
            root=work,
            datatype="surf",
            suffix="coords.nii.gz",
            space="corobl",
            desc="formesh",
            hemi="{hemi}",
            label="{label}",
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
            root=work,
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
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
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
            root=work,
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


# --- WIP rule for adjusting mesh in unfolded space


rule unfold_spring_model:
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
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldspring",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    script:
        "../scripts/unfold_spring_model.py"


# --- calculating metrics on native anatomical surfaces


rule calculate_surface_area:
    input:
        gii=bids(
            root=work,
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


rule calculate_gyrification:
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
            root=work,
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


rule calculate_curvature_from_surface:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=work,
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


rule calculate_thickness_from_surface:
    input:
        inner=bids(
            root=work,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        outer=bids(
            root=work,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        gii=bids(
            root=work,
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
            root=root,
            space="unfold",
            label="{label}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
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
        #"c3d {input.ref_nii} -scale 0 -shift 1 -pad 16x16x0vox 16x16x0vox 0 -o {output.ref_nii}"
        "c3d {input.ref_nii} -scale 0 -shift 1 -pad 64x64x0vox 64x64x0vox 0 -o {output.ref_nii}"


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


# TODO - if we add density wildcard here then it will make use of the template unfold standard densities
rule native_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs"""
    input:
        metric_gii=bids(
            root=work,
            datatype="surf",
            suffix="{metric}.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        inner_surf=bids(
            root=work,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        midthickness_surf=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
        outer_surf=bids(
            root=work,
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


# TODO - if we add density wildcard for inputs here, then it will make use of the template unfold standard densities,
# if no density wildcard, then will use the native mesh
#
# native metric, resampled to a standard density


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


rule unfoldreg:
    """performs 2D registration in unfolded space as in [reference paper]
    this is done in a shadow directory to get rid of the tmp files generated by ants."""
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
        "antsRegistrationSyNQuick.sh {params.antsparams} {params.fixed_args} {params.moving_args}  &> {log} && "
        "{params.cmd_copy_warps}"


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
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


def get_unfold_ref(wildcards):
    """function to return either unfoldreg or unfold ref mesh, depending on whether
    unfoldreg can be performed (based on atlas wildcards)"""

    if (
        wildcards.label in config["atlas_files"][config["atlas"]]["label_wildcards"]
        and config["no_unfolded_reg"] == False
    ):
        return bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfoldreg",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        )
    else:
        return bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        )


rule warp_unfold_native_to_unfoldreg:
    """ TODO: verify that this transformation is being correctly applied! """
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
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="native",
            to=config["atlas"],
            space="unfold",
            type_="surface",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
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


# --- resampling using the unfoldreg surface to standard densities (0p5mm, 1mm, 2mm, unfoldiso)
#
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
        None
    #        config["singularity"]["autotop"]
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
        None
    #        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-resample {input.native} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check"


rule resample_native_metric_to_std_density:
    input:
        native_metric=bids(
            root=work,
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
        None
    #        config["singularity"]["autotop"]
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
            root=work,
            datatype="surf",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container:
        None
    #        config["singularity"]["autotop"]
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
        ALSO - this has issues (with ribbon mapping) since the most inner and outer coords are missing.. 
             -- maybe nearest-vertex works better? 
"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_nii=bids(
            root=root,
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
        #use really large distance to ensure all voxels labelled


#        " -ribbon-constrained {params.inner_surf} {params.outer_surf}"
