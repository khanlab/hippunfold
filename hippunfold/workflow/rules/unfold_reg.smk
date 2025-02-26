# --- unfolded registration


# for unfold_reg, need to get data into 2d image - Jordan's previous code used metric-to-volume-mapping with unfoldiso inner and outer surfaces
# we can use same approach, just calculating the metrics on native, and using the unfold inner/outer as ribbon
unfoldreg_method = "SyN"  # choices: ["greedy","SyN"]

unfoldreg_padding = "64x64x0vox"


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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/c3d.yaml"
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input.ref_3d_nii} -slice z 50% -o {output.ref_2d_nii}"


rule space_unfold_vertices_evenly:
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
            **inputs.subj_wildcards,
        ),
    params:
        step_size=.1,
        max_iterations=10000,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    group:
        "subj"
    log:
        surf_gii=bids(
            root="logs",
            suffix="{surfname}.txt",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    script:
        "../scripts/space_unfold_vertices.py"


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
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_nii=bids(
            root=work,
            datatype="anat",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        interp="-nearest-vertex 0.3",
    output:
        metric_nii=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
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
            **inputs.subj_wildcards,
        ),
        metric=unfoldreg_metrics,
        **wildcards,
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
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="unfoldreg.txt",
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
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
    conda:
        "../envs/greedy.yaml"
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="unfoldreg_greedy.txt",
            atlas="{atlas}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
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
            type_="itk",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/neurovis.yaml"
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
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule warp_unfold_native_to_unfoldreg:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
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
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldreg",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shadow:
        "minimal"
    shell:
        "wb_command -volume-to-surface-mapping {input.warp} {input.surf_gii} warp.shape.gii -trilinear && "
        "wb_command -surface-coordinates-to-metric {input.surf_gii} coords.shape.gii && "
        "wb_command -metric-math 'COORDS + WARP' warpedcoords.shape.gii -var COORDS coords.shape.gii -var WARP warp.shape.gii && "
        "wb_command -surface-set-coordinates  {input.surf_gii} warpedcoords.shape.gii {output.surf_gii}"
