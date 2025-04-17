
rule gen_native_mesh:
    input:
        coords=lambda wildcards: bids(
            root=root,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="coords.nii.gz",
            desc=config["laminar_coords_method"],
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        nan_mask=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="nan",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        sink_mask=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="sink",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        src_mask=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="src",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        threshold=lambda wildcards: config["surf_thresholds"][wildcards.surfname],
        decimate_opts=0.75,
        hole_fill_radius=1.0,
    output:
        surf_gii=temp(
            temp(
                bids(
                    root=root,
                    datatype="surf",
                    suffix="{surfname,midthickness}.surf.gii",
                    den="native",
                    space="corobl",
                    desc="nostruct",
                    hemi="{hemi}",
                    label="{label}",
                    **inputs.subj_wildcards,
                )
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("pyvista")
    log:
        bids_log(
            "gen_native_mesh",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            desc="{surfname}",
        ),
    script:
        "../scripts/gen_isosurface.py"


rule update_native_mesh_structure:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="{space}",
            den="native",
            desc="nostruct",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname,midthickness|inner|outer}.surf.gii",
            space="{space,corobl|unfold}",
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
        "cp {input} {output} && wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule smooth_surface:
    """ slight smoothing of surface to improve curvature estimation """
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        smoothing_strength=0.8,
        smoothing_iterations=10,
    output:
        surf_gii=temp(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="corobl",
                den="native",
                desc="smoothed",
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
        "wb_command -surface-smoothing {input} {params.smoothing_strength} {params.smoothing_iterations} {output}"
