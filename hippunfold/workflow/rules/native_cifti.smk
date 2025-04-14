def get_inputs_cifti_metric_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_metric"] = (
            bids(
                root=root,
                datatype="metric",
                suffix="{metric}.shape.gii",
                den="native",
                hemi="L",
                label="{label}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_metric"] = (
            bids(
                root=root,
                datatype="metric",
                suffix="{metric}.shape.gii",
                den="native",
                hemi="R",
                label="{label}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dscalar_metric_cifti_native:
    input:
        unpack(get_inputs_cifti_metric_native),
    params:
        cmd=get_cmd_cifti_metric,
    output:
        cifti=bids(
            root=root,
            datatype="cifti",
            suffix="{metric}.dscalar.nii",
            den="native",
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
        "{params.cmd}"


def get_inputs_cifti_label_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_label"] = (
            bids(
                root=root,
                datatype="metric",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                den="native",
                hemi="L",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_label"] = (
            bids(
                root=root,
                datatype="metric",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                den="native",
                hemi="R",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dlabel_cifti_subfields_native:
    input:
        unpack(get_inputs_cifti_label_native),
    params:
        cmd=get_cmd_cifti_label,
    output:
        cifti=bids(
            root=root,
            datatype="cifti",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
            den="native",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "{params.cmd}"
