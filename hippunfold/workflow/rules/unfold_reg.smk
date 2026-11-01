
rule create_unfold_ref_slice:
    """Creates a single-slice 2D unfold ref at the midthickness z-level, with the same
    in-plane geometry as the 3D unfold ref. This is created directly rather than taking
    the 50% slice of the 3D ref, since that slice is offset by half a voxel from the
    midthickness plane, which made sampling of the 2D warps at the surface fail"""
    output:
        ref_2d_nii=temp(
            bids(
                root=root,
                datatype="warps",
                hemi="{hemi}",
                suffix="refvol.nii.gz",
                space="unfold",
                desc="slice",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        dims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["dims"][:2]
        )
        + "x1",
        voxdims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["voxdims"]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["origin"][:2]
        )
        + "x"
        + str(-1.0 * get_unfold_midthickness_z(wildcards.label)),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.label]["orient"],
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} {params.flip_per_hemi} -o {output.ref_2d_nii}"


rule native_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs"""
    input:
        metric_gii=bids(
            root=root,
            datatype="metric",
            suffix="{metric}.shape.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_nii=bids(
            root=root,
            datatype="warps",
            hemi="{hemi}",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        metric_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric,[0-9a-zA-Z]+}.nii.gz",
                space="unfold",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        interp="-nearest-vertex 10",
    shell:
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.midthickness_surf} {input.ref_nii} {output.metric_nii} "
        " {params.interp}"


rule atlas_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs.
    This rule is for the surface template"""
    input:
        ref_nii=bids(
            root=root,
            datatype="warps",
            hemi="{hemi}",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            suffix="{metric}.shape.gii",
        ),
        midthickness_surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="midthickness.surf.gii",
        ),
        inner_surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="inner.surf.gii",
        ),
        outer_surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="outer.surf.gii",
        ),
    output:
        metric_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold",
                hemi="{hemi}",
                label="{label}",
                den="{density}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    shell:
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.midthickness_surf} {input.ref_nii} {output.metric_nii} "
        " -ribbon-constrained {input.inner_surf} {input.outer_surf}"


rule slice_3d_to_2d_subject:
    """This is needed so ants will believe the data is truly 2d"""
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="{unfold}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        img=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="{unfold,unfold|unfoldLR}2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    params:
        clip_min=-5.0,
        clip_max=+5.0,
    script:
        "../scripts/slice_3d_to_2d.py"


rule slice_3d_to_2d_atlas:
    """This is needed so ants will believe the data is truly 2d"""
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    output:
        img=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                den="{density,[0-9k]+}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(1),
    params:
        clip_min=-5.0,
        clip_max=+5.0,
    script:
        "../scripts/slice_3d_to_2d.py"


def subject_unfoldreg_metric(metric, atlas):
    """Atlases written before the smoothing-based gyrification index store the
    legacy area ratio as gyrification, so the subject's legacy map is the one
    to register against them."""
    if (
        metric == "gyrification"
        and config["atlas_metadata"][atlas].get("gyrification_definition") != "sgi"
    ):
        return "legacygyrification"
    return metric


def get_fixed_images_unfoldreg(wildcards):
    unfoldreg_metrics = [
        subject_unfoldreg_metric(metric, wildcards.atlas)
        for metric in config["atlas_metadata"][wildcards.atlas]["metric_wildcards"]
    ]

    return expand(
        bids(
            root=root,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        metric=unfoldreg_metrics,
        **wildcards,
    )


def get_moving_images_unfoldreg(wildcards):
    unfoldreg_metrics = config["atlas_metadata"][wildcards.atlas]["metric_wildcards"]
    return expand(
        bids(
            root=root,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            den=config["unfoldreg_density"],
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
        metric=unfoldreg_metrics,
        **wildcards,
    )


def get_unfoldreg_syn_args(wildcards, input):
    """SyN stage for unfoldreg_ants, from unfoldreg_ants in the config"""
    params = config["unfoldreg_ants"]
    check_ants_levels(params)
    fixed, moving = input.fixed_images, input.moving_images
    weights = str(params["wt"]).split("x") if params.get("wt") else ["1"] * len(fixed)
    if len(weights) != len(fixed):
        raise ValueError(
            f"unfoldreg_ants weights {params['wt']} must have {len(fixed)} components"
        )
    metrics = " ".join(
        f"-m '{get_ants_syn_metric(params['sim'], f, m, w)}'"
        for f, m, w in zip(fixed, moving, weights)
    )
    return (
        f"-t '{get_ants_syn_transform(params)}' {metrics} "
        f"-c '[{params['q']},1e-9,10]' -f {params['f']} -s {params['s']}vox"
    )


rule unfoldreg_ants:
    """performs 2D registration in unfolded space as in [reference paper]
    this is done in a shadow directory to get rid of the tmp files generated by ants.

    Note: fixed and moving are swapped as compared to v1 unfoldreg.

    The SyN stage (unfoldreg_ants in the config) is the same as the one used by
    antsMultivariateTemplateConstruction2.sh to register the atlas inputs to the template
    (new_atlas_ants), with one metric per surface metric."""
    input:
        fixed_images=get_fixed_images_unfoldreg,
        moving_images=get_moving_images_unfoldreg,
    output:
        warp=temp(
            bids(
                root=root,
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
            )
        ),
        invwarp=temp(
            bids(
                root=root,
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
            )
        ),
    log:
        bids_log(
            "unfoldreg_ants",
            **inputs.subj_wildcards,
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
        ),
    group:
        "surf"
    shadow:
        "minimal"
    threads: 4
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(5),
    params:
        # A centre-of-mass init (-r [f,m,1], the antsRegistrationSyNQuick.sh default)
        # translates by the intensity-weighted centre of mass of the first metric
        # (curvature). Curvature is signed and roughly zero-mean,
        # so that centre is unstable and gives a spurious shift that SyN-only can't undo.
        # Subject and atlas maps share the same unfolded grid, so use the geometric
        # centre (-r [f,m,0]), which gives a zero initial translation.
        init_args=lambda wildcards, input: "-r '[{fixed},{moving},0]'".format(
            fixed=input.fixed_images[0], moving=input.moving_images[0]
        ),
        syn_args=get_unfoldreg_syn_args,
        cmd_copy_warps=lambda wildcards, output: " && ".join(
            [
                f"cp tmp1Warp.nii.gz {output.warp}",
                f"cp tmp1InverseWarp.nii.gz {output.invwarp}",
            ]
        ),
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsRegistration -d 2 --float 1 --verbose 1 -u 0 -w '[0.01,0.99]' -z 1 -o tmp "
        "{params.init_args} {params.syn_args} &> {log} && "
        "{params.cmd_copy_warps}"


rule reset_header_2d_warp_unfoldreg:
    """adjusts header to match the original data
    (since this seems to get garbled in z by ants)"""
    input:
        nii=bids(
            root=root,
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
        ref_nii=bids(
            root=root,
            datatype="warps",
            hemi="{hemi}",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=temp(
            bids(
                root=root,
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
            )
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(1),
    script:
        "../scripts/set_metric_nii_header.py"


rule warp_unfold_native_to_unfoldreg:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_=config["atlas"],
            to="native",
            space="unfold",
            type_="itk",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        surf_gii=temp(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="unfoldreg",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "surf"
    shadow:
        "minimal"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=get_cmd_warp_surface_2d_warp,
    shell:
        "{params.cmd}"
