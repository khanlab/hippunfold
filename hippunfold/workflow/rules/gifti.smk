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
        structure_type=lambda wildcards: get_structure(
            wildcards.hemi, wildcards.autotop
        ),
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
            **inputs.subj_wildcards,
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
        structure_type=lambda wildcards: get_structure(
            wildcards.hemi, wildcards.autotop
        ),
    output:
        coords_gii=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="coords.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
        ),
        ref_nii=bids(
            root=work,
            datatype="warps",
            space="unfold",
            label="{autotop}",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
        ),
    log:
        bids(
            root="logs",
            den="{density}",
            suffix="{surfname}.txt",
            desc="constrainbbox",
            hemi="{hemi}",
            label="{autotop}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/constrain_surf_to_bbox.py"


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
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=work,
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
            label="{autotop,hipp|dentate}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
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
                label="{autotop}",
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
                label="{autotop}",
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
            label="{autotop}",
            **inputs.subj_wildcards,
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
                **inputs.subj_wildcards,
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
            den="{density}",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
            space="{space}",
            label="hipp",
            **inputs.subj_wildcards,
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
        structure = get_structure(wildcards.hemi, wildcards.label)
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
                den="{density}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
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
            space=["{space}", "unfold", "unfoldreg"],
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
                den="{density}",
                suffix="subfields.dlabel.nii",
                atlas="{atlas}",
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
                label="{autotop}",
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
            label="{autotop}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
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
                label="{autotop}",
                **inputs.subj_wildcards,
            ),
            autotop=config["autotop_labels"],
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
    group:
        "subj"
    shell:
        "{params.cmd}"
