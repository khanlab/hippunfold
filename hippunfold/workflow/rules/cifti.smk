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
                datatype="metric",
                den="{density}",
                suffix="{metric}.shape.gii",
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
                den="{density}",
                suffix="{metric}.shape.gii",
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
            datatype="cifti",
            den="{density}",
            suffix="{metric}.dscalar.nii",
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
                datatype="metric",
                den="{density}",
                atlas="{atlas}",
                suffix="subfields.label.gii",
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
                den="{density}",
                atlas="{atlas}",
                suffix="subfields.label.gii",
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
            datatype="cifti",
            den="{density}",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
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


rule create_spec_file:
    input:
        lambda wildcards: get_inputs_spec_file(wildcards.label, wildcards.density),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=temp(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="surfaces.spec",
                hemi="{hemi,L|R}",
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
        "{params.cmds}"


def get_cmd_merge_spec(wildcards, input, output):
    if len(input.spec_files) == 1:
        return f"cp {input} {output}"
    else:
        return f"wb_command -spec-file-merge {input.spec_files} {output}"


rule merge_lr_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="{surfdir}",
                den="{density}",
                suffix="surfaces.spec",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=temp(
            bids(
                root=root,
                datatype="{surfdir}",
                den="{density}",
                suffix="surfaces.spec",
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
        "{params.cmd}"


rule merge_hipp_dentate_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="{surfdir}",
                den="{density}",
                suffix="surfaces.spec",
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
            datatype="{surfdir}",
            den="{density}",
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


def get_inputs_to_remove(wildcards):
    files = []
    for label in config["autotop_labels"]:
        for spec_input in get_inputs_spec_file(label, density=config["unused_density"]):
            files.extend(
                inputs[config["modality"]].expand(
                    spec_input,
                    label=label,
                    **wildcards,
                    **expand_hemi(),
                )
            )
    files.extend(
        expand(
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
            surfname=["midthickness", "inner", "outer"],
            space=["corobl", "unfold"],
            label=config["autotop_labels"],
            density=config["unused_density"],
            **wildcards,
            **expand_hemi(),
            allow_missing=True,
        )
    )
    return files


rule remove_extra_spec_inputs:
    """ after all specs are created, remove files with unused density """
    input:
        expand(
            bids(
                root=root,
                datatype="{surfdir}",
                den="{density}",
                suffix="surfaces.spec",
                **inputs.subj_wildcards,
            ),
            density=config["output_density"],
            space=ref_spaces,
            allow_missing=True,
        ),
    params:
        to_remove=get_inputs_to_remove,
    output:
        temp(
            bids(
                root=root,
                datatype="{surfdir}",
                suffix="removeunused.touch",
                **inputs.subj_wildcards,
            )
        ),
    shell:
        "rm -f {params.to_remove} && touch {output}"
