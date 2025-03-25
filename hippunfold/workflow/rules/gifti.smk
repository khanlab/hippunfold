# lookup mappings for structure and surface type:
def get_structure(hemi, label):
    if label == "hipp":
        if hemi == "L":
            return "HIPPOCAMPUS_LEFT"
        elif hemi == "R":
            return "HIPPOCAMPUS_RIGHT"
    elif label == "dentate":
        if hemi == "L":
            return "HIPPOCAMPUS_DENTATE_LEFT"
        elif hemi == "R":
            return "HIPPOCAMPUS_DENTATE_RIGHT"


surf_to_secondary_type = {
    "midthickness": "MIDTHICKNESS",
    "inner": "INNER",
    "outer": "OUTER",
}


# warp from corobl to native
rule affine_gii_to_native:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="ras",
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="{native_modality,T1w|T2w}",
            hemi="{hemi}",
            label="{label,hipp|dentate}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("ants")
    group:
        "subj"
    shell:
        "wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}"


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
                label="{label}",
                **inputs.subj_wildcards,
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
                label="{label}",
                **inputs.subj_wildcards,
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


def get_inputs_cifti_label(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_label"] = (
            bids(
                root=root,
                datatype="surf",
                den="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="L",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_label"] = (
            bids(
                root=root,
                datatype="surf",
                den="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="R",
                label="hipp",
                **inputs.subj_wildcards,
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
            den="{atlas}",
            suffix="subfields.dlabel.nii",
            space="{space}",
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


def get_cmd_spec_file(wildcards, input, output):
    specfile = output.spec_file
    if "hemi" in wildcards._names:
        structure = get_structure(wildcards.hemi, wildcards.label)
    else:
        structure = "INVALID"
    cmds = list()
    for infile in input:
        cmds.append(
            " ".join(["wb_command", "-add-to-spec-file", specfile, structure, infile])
        )
    return " && ".join(cmds)


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
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        subfields=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                den="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
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
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", "unfoldreg"],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        cifti_labels=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                den="{atlas}",
                suffix="subfields.dlabel.nii",
                space="{space}",
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
            datatype="surf",
            den="{density}",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
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
                **inputs.subj_wildcards,
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
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", "unfoldreg"],
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
            datatype="surf",
            den="{density}",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
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


def get_cmd_merge_spec(wildcards, input, output):
    if len(input.spec_files) == 1:
        return f"cp {input} {output}"
    else:
        return f"wb_command -spec-file-merge {input.spec_files} {output}"


rule merge_lr_spec_file_native:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="surfaces.spec",
                hemi="{hemi}",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            den="{density}",
            space="{space}",
            suffix="surfaces.spec",
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


rule merge_hipp_dentate_spec_file_native:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="surfaces.spec",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=config["autotop_labels"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            den="{density}",
            space="{space}",
            suffix="surfaces.spec",
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
