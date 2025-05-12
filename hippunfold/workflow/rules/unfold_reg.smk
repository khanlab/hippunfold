
# unfold ref nifti
rule create_unfold_ref:
    params:
        dims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["dims"]
        ),
        voxdims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["voxdims"]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["origin"]
        ),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.label]["orient"],
    output:
        nii=temp(
            bids(
                root=root,
                space="unfold",
                label="{label}",
                datatype="warps",
                suffix="refvol.nii.gz",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}"



rule extract_unfold_ref_slice:
    """This gets the central-most slice of the unfold volume, for obtaining a 2D slice"""
    input:
        ref_3d_nii=bids(
            root=root,
            datatype="warps",
            suffix="refvol.nii.gz",
            space="unfold",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        ref_2d_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="refvol.nii.gz",
                space="unfold",
                desc="slice",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input.ref_3d_nii} -slice z 50% -o {output.ref_2d_nii}"


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
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        interp="-nearest-vertex 10",
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
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.midthickness_surf} {input.ref_nii} {output.metric_nii} "
        " {params.interp}"


rule atlas_metric_to_unfold_nii:
    """converts metric .gii files to .nii for use in ANTs. 
        This rule is for the surface template"""
    input:
        ref_nii=bids(
            root=root,
            datatype="anat",
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
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
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
            space="unfold",
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
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        conda_env("neurovis")
    group:
        "subj"
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
    conda:
        conda_env("neurovis")
    group:
        "subj"
    script:
        "../scripts/slice_3d_to_2d.py"


def get_fixed_images_unfoldreg(wildcards):
    unfoldreg_metrics = config["atlas_metadata"][wildcards.atlas]["metric_wildcards"]

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
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    log:
        bids_log(
            "unfoldreg_antsquick",
            **inputs.subj_wildcards,
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
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


rule reset_header_2d_warp_unfoldreg:
    """ adjusts header to match the original data
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
            suffix="refvol.nii.gz",
            space="unfold",
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
    conda:
        conda_env("neurovis")
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
    params:
        cmd=get_cmd_warp_surface_2d_warp,
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
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shadow:
        "minimal"
    shell:
        "{params.cmd}"
