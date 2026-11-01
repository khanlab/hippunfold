
rule calculate_surface_area:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="{space}",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="surfarea{space}.shape.gii",
                den="native",
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
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule metric_smoothing:
    input:
        surface=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        metric=bids(
            root=root,
            datatype="metric",
            suffix="{metric}.shape.gii",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        metric=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="{metric}.shape.gii",
                den="{density}",
                desc="fwhm{fwhm}mm",
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
        fwhm=lambda wildcards: str(wildcards.fwhm).replace("p", "."),
    shell:
        "wb_command -metric-smoothing {input.surface} {input.metric} {params.fwhm} {output.metric} -fwhm"


rule calculate_legacy_gyrification:
    """Ratio of native to unfolded (spring model) vertex area, each smoothed
    at 1 mm FWHM. Superseded by the smoothing-based gyrification index; kept
    because existing unfoldreg atlases were built with it."""
    input:
        native_surfarea=bids(
            root=root,
            datatype="metric",
            suffix="surfareacorobl.shape.gii",
            den="native",
            desc="fwhm1mm",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        unfold_surfarea=bids(
            root=root,
            datatype="metric",
            suffix="surfareaunfoldspringmodelsmooth.shape.gii",
            den="native",
            desc="fwhm1mm",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="legacygyrification.shape.gii",
                desc="unnorm",
                den="native",
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
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea}"


rule calculate_gyrification:
    """Smoothing-based gyrification index (sGI): the midthickness area over
    the area of a heat-smoothed copy of itself (boundary pinned), pooled on
    the surface. Folds narrower than about the scale count; see
    lib/gyrification.py."""
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="gyrification.shape.gii",
                desc="unnorm",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    log:
        bids_log(
            "calculate_gyrification",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(5),
    params:
        scale=config["gyrification_scale"],
        sigma=config["gyrification_pooling_sigma"],
    script:
        "../scripts/calculate_gyrification.py"


rule calculate_curvature:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="native",
            desc="smoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="curvature.shape.gii",
                den="native",
                desc="unnorm",
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
    shell:
        "wb_command -surface-curvature {input} -mean {output}"


rule calculate_thickness:
    input:
        inner=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="thickness.shape.gii",
                den="native",
                desc="unnorm",
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
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"


rule soft_tanh_normalization:
    """Normalize the surface shape metrics (thick, curv, gyr, legacy gyr) to smoothly
    reject outliers that can cause issues with registration"""
    input:
        gii=bids(
            root=root,
            datatype="metric",
            suffix="{metric}.shape.gii",
            den="native",
            desc="unnorm",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="{metric,thickness|curvature|gyrification|legacygyrification}.shape.gii",
                den="native",
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
        norm_limit=5.0,
    shell:
        "wb_command -metric-math '{params.norm_limit}*tanh(X/{params.norm_limit})' {output.gii} -var X {input.gii}"
