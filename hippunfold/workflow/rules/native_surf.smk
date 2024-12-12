surf_thresholds={'inner': 0.05, 'outer':0.95, 'midthickness':0.5}
print(bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="unfold",
            desc="{desc}nostruct",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards
        ))

print(bids(
            root=root,
            datatype="surf_",
            fromdensity="{density}",
            suffix="subfields.label.gii",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ))

rule gen_native_hipp_mesh:
    input:
        nii=bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="hipp",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    params:
        threshold=lambda wildcards: surf_thresholds[wildcards.surfname],
        decimate_percent=0
    output:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}nostruct",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    #container (will need pyvista dependency)
    script:
        "../scripts/gen_isosurface.py"

rule update_native_hipp_mesh_structure:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}nostruct",
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
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}",
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
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="{desc}",
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
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="{surfname}.surf.gii",
            space="unfold",
            desc="{desc}",
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
            root=root,
            datatype="surf",
            den="{density}",
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
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=bids(
            root=root,
            datatype="surf_", #temporarily, to keep things separate for development
            suffix="midthickness.surf.gii",
            space="unfold",
            desc="{desc}",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf_",
            fromdensity="{density}",
            desc="{desc}",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **inputs.subj_wildcards
        ),
    container: None
#        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        'wb_command -label-resample {input.label_gii} {input.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check'

