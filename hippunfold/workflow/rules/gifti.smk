# lookup tables for structure:
hemi_to_structure = {"L": "CORTEX_LEFT", "Lflip": "CORTEX_LEFT", "R": "CORTEX_RIGHT"}
surf_to_secondary_type = {
    "midthickness": "MIDTHICKNESS",
    "inner": "PIAL",
    "outer": "GRAY_WHITE",
}


rule cp_template_to_unfold:
    """cp template unfold surf to subject"""
    input:
        gii=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{autotop}",
            "tpl-avg_space-unfold_den-{density}_{surfname}.surf.gii",
        ),
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="FLAT",
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "cp {input.gii} {output.gii} && "
        "wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule calc_unfold_template_coords:
    """ Creates a coords.shape.gii from the unfolded template """
    input:
        midthickness_gii=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
    params:
        coords_xyz="coords-XYZ.shape.gii",
        coord_AP="coord-AP.shape.gii",
        coord_PD="coord-PD.shape.gii",
        coord_IO="coord-IO.shape.gii",
        origin=lambda wildcards: config["{autotop}"]["origin"],
        extent=lambda wildcards: config["{autotop}"]["extent"],
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
    output:
        coords_gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="coords.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"  #this is required to use the temporary files defined as params
    group:
        "subj"
    shell:
        "wb_command -surface-coordinates-to-metric {input.midthickness_gii} {params.coords_xyz} && "
        "wb_command -file-information {params.coords_xyz} && "
        "wb_command -metric-math '({params.origin[0]}-X)/{params.extent[0]}' {params.coord_AP} -var X {params.coords_xyz} -column 1 && "
        "wb_command -file-information {params.coord_AP} && "
        "wb_command -metric-math '({params.origin[1]}+Y)/{params.extent[1]}' {params.coord_PD} -var Y {params.coords_xyz} -column 2 && "
        "wb_command -file-information {params.coord_PD} && "
        "wb_command -metric-math '({params.origin[2]}+Z)/{params.extent[2]}' {params.coord_IO} -var Z {params.coords_xyz} -column 3 && "
        "wb_command -file-information {params.coord_IO} && "
        "wb_command -metric-merge {output.coords_gii}  -metric {params.coord_AP} -metric {params.coord_PD} -metric {params.coord_IO} && "
        "wb_command -set-structure {output.coords_gii} {params.structure_type}"


# subj unfolded surf might have a few vertices outside the bounding box.. this constrains all the vertices to the warp bounding box
rule constrain_surf_to_bbox:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        ref_nii=bids(
            root=root,
            datatype="warps",
            space="unfold",
            label="{autotop}",
            suffix="refvol.nii.gz",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    log:
        bids(
            root="logs",
            den="{density}",
            suffix="{surfname}.txt",
            desc="constrainbbox",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/constrain_surf_to_bbox.py"


### A) we will run the following rules once for unfoldreg-none and again later (with no unfoldreg wildcard, as in previous version)


# warp from subj unfolded to corobl
rule warp_gii_unfold2corobl1:
    input:
        warp=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            label="hipp",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="surface"
        ),
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nonancorrect",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && "
        "wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


# previous rule seems to be where nan vertices emerge, so we'll correct them here immediately after
rule correct_nan_vertices1:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nonancorrect",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/fillnanvertices.py"


# morphological features, calculated in native space:
rule calculate_surface_area1:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="surfarea.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule calculate_gyrification1:
    """new gyrification is ratio of nativearea to unfoldarea (e.g. surface scaling or distortion factor.
    this should be proportional by a constant, to the earlier gyrification on 32k surfaces."""
    input:
        native_surfarea=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="surfarea.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        unfold_surfarea=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_hipp",
            "tpl-avg_space-unfold_den-{density}_surfarea.shape.gii",
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="gyrification.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    log:
        bids(
            root="logs",
            den="{density}",
            suffix="calcgyrification.txt",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea} &> {log}"


rule calculate_curvature_from_surface1:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            desc="unnorm",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-curvature {input} -mean {output}"


rule normalize_curvature1:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            desc="unnorm",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/normalize_tanh.py"


rule calculate_thickness_from_surface1:
    input:
        inner=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="inner.surf.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        outer=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="outer.surf.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="thickness.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"


### B) up to here everything is calculated in unfoldreg-none. Now we will align to the reference atlas based on thickness, curvature, and gyrification


rule metric_to_nii:
    """converts metric .gii files to .nii for use in ANTs"""
    input:
        metric_gii=bids(
            root=work,
            datatype="surf",
            den="unfoldiso",
            suffix="{metric}.shape.gii",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        unfoldedsurf=bids(
            root=work,
            datatype="surf",
            den="unfoldiso",
            suffix="inner.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        atlas_dir=Path(download_dir) / "atlas" / "multihist7",
    params:
        interp="-nearest-vertex 1",
        refflatnii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"]["multihist7"]["label_nii"].format(**wildcards),
    output:
        metric_nii=bids(
            root=work,
            datatype="anat",
            den="unfoldiso",
            suffix="{metric}.nii.gz",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -metric-to-volume-mapping {input.metric_gii} {input.unfoldedsurf} {params.refflatnii} {output.metric_nii} {params.interp}"


rule unfolded_registration:
    """performs 2D registration in unfolded space as in [reference paper]
    this is done in a shadow directory to get rid of the tmp files generated by ants."""
    input:
        thickness=bids(
            root=work,
            datatype="anat",
            den="unfoldiso",
            suffix="thickness.nii.gz",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        curvature=bids(
            root=work,
            datatype="anat",
            den="unfoldiso",
            suffix="curvature.nii.gz",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        gyrification=bids(
            root=work,
            datatype="anat",
            den="unfoldiso",
            suffix="gyrification.nii.gz",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
    params:
        antsparams="-d 2 -t so",
        outsuffix="tmp",
        warpfn="tmp1Warp.nii.gz",
        invwarpfn="tmp1InverseWarp.nii.gz",
        refthickness=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["thickness"].format(**wildcards),
        refcurvature=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["curvature"].format(**wildcards),
        refgyrification=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["gyrification"].format(**wildcards),
    output:
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="{atlas}",
            to="subject",
            space="unfold",
            type_="itk",
            hemi="{hemi}"
        ),
        invwarp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="subject",
            to="{atlas}",
            space="unfold",
            type_="itk",
            hemi="{hemi}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    log:
        bids(
            root="logs",
            den="unfoldiso",
            suffix="unfoldedRegistration.txt",
            space="unfold",
            unfoldreg="{atlas}",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    shadow:
        "minimal"
    threads: 16
    resources:
        mem_mb=16000,
        time=10,
    shell:
        "antsRegistrationSyNQuick.sh {params.antsparams} -f {params.refthickness} -f {params.refcurvature} -f {params.refgyrification} -m {input.thickness} -m {input.curvature} -m {input.gyrification} -o {params.outsuffix} &> {log} && "
        "cp {params.warpfn} {output.warp} && "
        "cp {params.invwarpfn} {output.invwarp}"


# warp from subj unfolded to unfoldedaligned
rule warp_gii_unfoldreg:
    input:
        invwarp=expand(
            bids(
                root=work,
                **config["subj_wildcards"],
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="SyN",
                from_="subject",
                to="{atlas}",
                space="unfold",
                type_="itk",
                hemi="{hemi}"
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/warp_flatsurf.py"


rule dentate_skip_unfoldreg:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="dentate",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            hemi="{hemi}",
            label="dentate",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


### C) We will now repeat A), then continue below


def skip_unfoldreg_option(wildcards):
    if config["no_unfolded_reg"]:
        gii = bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            unfoldreg="none",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        )
    else:
        gii = bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="constrainbbox",
            space="unfold",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        )
    return gii


# warp from subj unfolded to corobl
rule warp_gii_unfold2corobl2:
    input:
        warp=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="surface"
        ),
        gii=skip_unfoldreg_option,
    params:
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nonancorrect",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && "
        "wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


# previous rule seems to be where nan vertices emerge, so we'll correct them here immediately after
rule correct_nan_vertices2:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nonancorrect",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/fillnanvertices.py"


# needed for if native_modality is corobl
rule cp_corobl_root:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


# warp from corobl to native
rule affine_gii_to_native:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
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
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="{native_modality}",
            hemi="{hemi}",
            label="{autotop,hipp|dentate}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}"


# morphological features, calculated in native space:
rule calculate_surface_area2:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="surfarea.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule calculate_gyrification2:
    """new gyrification is ratio of nativearea to unfoldarea (e.g. surface scaling or distortion factor.
    this should be proportional by a constant, to the earlier gyrification on 32k surfaces."""
    input:
        native_surfarea=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="surfarea.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        unfold_surfarea=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{autotop}",
            "tpl-avg_space-unfold_den-{density}_surfarea.shape.gii",
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="gyrification.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    log:
        bids(
            root="logs",
            den="{density}",
            suffix="calcgyrification.txt",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea} &> {log}"


rule calculate_curvature_from_surface2:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="{space}",
            hemi="{hemi}",
            desc="unnorm",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-curvature {input} -mean {output}"


rule normalize_curvature2:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="{space}",
            hemi="{hemi}",
            desc="unnorm",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="curvature.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/normalize_tanh.py"


rule calculate_thickness_from_surface2:
    input:
        inner=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="inner.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        outer=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="outer.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="thickness.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"


### D) now continue as in previous version


rule resample_atlas_to_refvol:
    """this is just done in case the atlas has a different unfolded config than the current run"""
    input:
        refvol=bids(
            root=root,
            space="unfold",
            label="hipp",
            datatype="warps",
            suffix="refvol.nii.gz",
            **config["subj_wildcards"]
        ),
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
    params:
        atlas=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_nii"].format(**wildcards),
    output:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    log:
        bids(
            root="logs",
            suffix="resamplesubfieldrefvol",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {params.atlas} -r {input.refvol} -o {output.label_nii} -v &> {log}"


rule nii_to_label_gii:
    input:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
        surf=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_hipp",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
    params:
        label_list=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_list"].format(**wildcards),
        structure_type=lambda wildcards: hemi_to_structure[wildcards.hemi],
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="subfields.label.gii",
            space="{space}",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"
    shell:
        "wb_command -volume-to-surface-mapping {input.label_nii} {input.surf} temp.shape.gii -enclosing  && "
        "wb_command -metric-label-import temp.shape.gii {params.label_list} {output.label_gii} && "
        "wb_command -set-structure {output.label_gii} {params.structure_type}"


def get_cmd_cifti_metric(wildcards, input, output):
    cmd = f"wb_command  -cifti-create-dense-scalar {output}"
    if "L" in config["hemi"]:
        cmd = cmd + f" -left-metric {input.left_metric}"
    if "R" in config["hemi"]:
        cmd = cmd + f" -right-metric {input.right_metric}"
    return cmd


def get_inputs_cifti_metric(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_metric"] = (
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="L",
                label="{autotop}",
                **config["subj_wildcards"],
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_metric"] = (
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="R",
                label="{autotop}",
                **config["subj_wildcards"],
            ).format(**wildcards),
        )
    return files


rule create_dscalar_metric_cifti:
    input:
        unpack(get_inputs_cifti_metric),
    params:
        cmd=get_cmd_cifti_metric,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{metric}.dscalar.nii",
            space="{space}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


def get_inputs_cifti_label(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_label"] = (
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="L",
                label="hipp",
                **config["subj_wildcards"],
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_label"] = (
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="R",
                label="hipp",
                **config["subj_wildcards"],
            ).format(**wildcards),
        )
    return files


def get_cmd_cifti_label(wildcards, input, output):
    cmd = f"wb_command  -cifti-create-label {output}"
    if "L" in config["hemi"]:
        cmd = cmd + f" -left-label {input.left_label}"
    if "R" in config["hemi"]:
        cmd = cmd + f" -right-label {input.right_label}"
    return cmd


rule create_dlabel_cifti_subfields:
    input:
        unpack(get_inputs_cifti_label),
    params:
        cmd=get_cmd_cifti_label,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            den="{density}",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
            space="{space}",
            label="hipp",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


def get_cmd_spec_file(wildcards, input, output):
    specfile = output.spec_file
    if "hemi" in wildcards._names:
        structure = hemi_to_structure[wildcards.hemi]
    else:
        structure = "INVALID"
    cmds = list()
    for infile in input:
        cmds.append(
            " ".join(["wb_command", "-add-to-spec-file", specfile, structure, infile])
        )
    return " && ".join(cmds)


def get_cifti_metric_types(label):
    types_list = config["cifti_metric_types"][label]
    if config["generate_myelin_map"]:
        types_list.append("myelin.dscalar")
    return types_list


def get_gifti_metric_types(label):
    types_list = config["gifti_metric_types"][label]
    if config["generate_myelin_map"]:
        types_list.append("myelin.shape")
    return types_list


rule cp_unfolded_noconstrain:
    input:
        gii=skip_unfoldreg_option,
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


rule create_spec_file_hipp:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **config["subj_wildcards"]
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        subfields=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **config["subj_wildcards"]
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
        surfs=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **config["subj_wildcards"]
            ),
            surfname=["midthickness"],
            space=["{space}", "unfold"],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **config["subj_wildcards"]
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        cifti_labels=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="subfields.dlabel.nii",
                atlas="{atlas}",
                space="{space}",
                label="{label}",
                **config["subj_wildcards"]
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
            den="{density}",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,hipp}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmds}"


rule create_spec_file_dentate:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **config["subj_wildcards"]
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        surfs=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **config["subj_wildcards"]
            ),
            surfname=["midthickness"],
            space=["{space}", "unfold"],
            allow_missing=True,
        ),
        cifti=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **config["subj_wildcards"]
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
            den="{density}",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,dentate}",
            **config["subj_wildcards"]
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
                den="{density}",
                suffix="surfaces.spec",
                hemi="{hemi}",
                space="{space}",
                label="{autotop}",
                **config["subj_wildcards"]
            ),
            hemi=["L", "R"],
            allow_missing=True,
        ),
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            den="{density}",
            space="{space}",
            suffix="surfaces.spec",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -spec-file-merge {input.spec_files} {output}"
