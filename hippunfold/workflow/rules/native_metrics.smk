
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
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
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
    params:
        fwhm=lambda wildcards: str(wildcards.fwhm).replace("p", "."),
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
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "wb_command -metric-smoothing {input.surface} {input.metric} {params.fwhm} {output.metric} -fwhm"


rule calculate_gyrification:
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
        gii=bids(
            root=root,
            datatype="metric",
            suffix="gyrification.shape.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea}"


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
        gii=bids(
            root=root,
            datatype="metric",
            suffix="curvature.shape.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
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
        gii=bids(
            root=root,
            datatype="metric",
            suffix="thickness.shape.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"
