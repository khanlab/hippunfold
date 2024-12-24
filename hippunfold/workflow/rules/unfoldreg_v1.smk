
### A) we will run the following rules once for unfoldreg-none and again later (with no unfoldreg wildcard, as in previous version)


# warp from subj unfolded to corobl
rule warp_gii_unfold2corobl1:
    input:
        warp=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards
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
            desc="nobadcorrect",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && "
        "wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


# previous rule seems to be where bad vertices emerge, so we'll correct them here immediately after
rule correct_bad_vertices1:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nobadcorrect",
            space="corobl",
            unfoldreg="none",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    params:
        dist=lambda wildcards: config["outlier_opts"]["outlierSmoothDist"][
            wildcards.density
        ],
        threshold=config["outlier_opts"]["vertexOutlierThreshold"],
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
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/fillbadvertices.py"


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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards
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
                **inputs.subj_wildcards,
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
        )
    return gii


# warp from subj unfolded to corobl
rule warp_gii_unfold2corobl2:
    input:
        warp=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
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
            desc="nobadcorrect",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && "
        "wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


# previous rule seems to be where bad vertices emerge, so we'll correct them here immediately after
rule correct_bad_vertices2:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            desc="nobadcorrect",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
    params:
        dist=lambda wildcards: config["outlier_opts"]["outlierSmoothDist"][
            wildcards.density
        ],
        threshold=config["outlier_opts"]["vertexOutlierThreshold"],
    output:
        gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/fillbadvertices.py"


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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
        ),
    log:
        bids(
            root="logs",
            den="{density}",
            suffix="calcgyrification.txt",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
        ),
        outer=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="outer.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"

### D) now continue as in previous version



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
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"



