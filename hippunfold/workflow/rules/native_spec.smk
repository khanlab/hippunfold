

rule create_spec_file_hipp_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="metric",
                suffix="{metric}.gii",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        subfields=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="metric",
                suffix="subfields.label.gii",
                den="native",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surfnative",
                suffix="{surfname}.surf.gii",
                space="{space}",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="cifti",
                suffix="{cifti}.nii",
                den="native",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        cifti_labels=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="cifti",
                suffix="subfields.dlabel.nii",
                atlas="{atlas}",
                den="native",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surfnative",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            den="native",
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "{params.cmds}"


rule create_spec_file_dentate_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="metric",
                suffix="{metric}.gii",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surfnative",
                suffix="{surfname}.surf.gii",
                space="{space}",
                den="native",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="cifti",
                suffix="{cifti}.nii",
                den="native",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surfnative",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            den="native",
            label="{label,dentate}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "{params.cmds}"
