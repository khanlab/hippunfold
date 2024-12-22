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


# if we aren't doing unfolded reg, then we have these rules override correct_bad_vertices2
ruleorder: resample_native_surf > correct_bad_vertices2
ruleorder: warp_midthickness_to_inout > correct_bad_vertices2


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


rule resample_subfields_to_native_surf:
    input:
        label_gii=bids(
            root=work,
            datatype="surf",
            den=resample_from_density,
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
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii".format(
                density=resample_from_density
            ),
        ),
        native_unfold=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
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
        "wb_command -label-resample {input.label_gii} {input.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


rule resample_native_surf:
    input:
        native_corobl=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
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
            "tpl-avg_space-unfold_den-{density}_{surf_name}.surf.gii",
        ),
        native_unfold=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        native_corobl_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="corobl",
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
        "wb_command -surface-resample {input.native_corobl} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_corobl_resampled} -bypass-sphere-check"


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


rule convert_warp_from_itk_to_world:
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
            den="{density}",
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
            den="{density}",
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
