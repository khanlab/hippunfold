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


unfoldreg_metrics = [ 'thickness', 'curvature', 'gyrification']


# if we aren't doing unfolded reg, then we have these rules override correct_bad_vertices2
ruleorder: resample_native_surf >  cp_template_to_unfold
#ruleorder: cp_template_to_unfold > resample_native_surf
ruleorder: resample_native_surf > correct_bad_vertices2 
ruleorder: warp_midthickness_to_inout > correct_bad_vertices2
ruleorder: unfold_spring_model > warp_unfold_native_to_unfoldreg

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
            "tpl-avg_space-unfold_den-{density}_{surf_name}.surf.gii",
        ),
        native_unfold=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    output:
        native_corobl_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
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
    shell:
        "greedy -d 3  -rf {input.fixed} -rm {input.moving} {output.warped}  -r {input.warp} "

#TODO: rename the warps here according to existing custom (type=itk -> type=ras)
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
    script: '../scripts/unfold_spring_model.py'



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




#for unfold_reg, need to get data into 2d image - Jordan's previous code used metric-to-volume-mapping with unfoldiso inner and outer surfaces
# we can use same approach, just calculating the metrics on native, and using the unfold inner/outer as ribbon

rule pad_unfold_ref:
    """pads the unfolded ref in XY (to improve registration by ensuring zero-displacement at boundaries)
        and to deal with input vertices that are just slightly outside the original bounding box"""
    input:
        atlas_dir=Path(download_dir) / "atlas" / "multihist7",
    params:
        ref_nii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"]["multihist7"]["label_nii"].format(**wildcards),
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
        #"c3d {params.ref_nii} -scale 0 -shift 1 -pad 16x16x0vox 16x16x0vox 0 -o {output.ref_nii}"
        "c3d {params.ref_nii} -scale 0 -shift 1 -pad 64x64x0vox 64x64x0vox 0 -o {output.ref_nii}"


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
        unfoldedsurf=bids(
            root=work,
            datatype="surf",
            suffix="inner.surf.gii", #inner is because the multihist 2d slices are there, but ideally this should be midthickness
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
            desc="padded",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ),
    params:
        interp="-nearest-vertex 1",
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
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.unfoldedsurf} {input.ref_nii} {output.metric_nii} {params.interp}"


print(bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="native",
            to_="{atlas}",
            space="unfold",
            type_="itk",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ))
    
rule unfoldreg:
    """performs 2D registration in unfolded space as in [reference paper]
    this is done in a shadow directory to get rid of the tmp files generated by ants."""
    input:
        fixed_images=expand(bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
            ),
            metric=unfoldreg_metrics,allow_missing=True),
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
    params:
        antsparams="-d 2 -t so -o tmp",
        fixed_args=lambda wildcards, input:
            ' '.join(['-f {img}'.format(img=img) for img in  input.fixed_images]),
        moving_args=lambda wildcards, input:
            ' '.join(['-m {img}'.format(img=str( Path(input.atlas_dir) / config["atlas_files"][wildcards.atlas][metric].format(**wildcards)))
                for metric in unfoldreg_metrics ] ),
        cmd_copy_warps=lambda wildcards, output:
            ' && '.join([f"cp tmp1Warp.nii.gz {output.warp}",
                            f"cp tmp1InverseWarp.nii.gz {output.invwarp}"])
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
        #"antsRegistrationSyNQuick.sh {params.antsparams} {params.fixed_args} {params.moving_args}  &> {log} && "
        "antsRegistrationSyNQuick.sh {params.antsparams} {params.fixed_args} {params.moving_args}  && "
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
    script: '../scripts/convert_warp_2d_to_3d.py'


rule convert_unfoldreg_warp_from_itk_to_world:
    input:
        warp=bids(
            root=work,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="{desc}",
            from_="{to}", #swapped around since surf needs opposite direction
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




rule warp_unfold_native_to_unfoldreg:
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
            to="{atlas}",
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
            space="unfold{atlas}",
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


       
