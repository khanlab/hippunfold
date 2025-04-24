from snakebids.paths import bids_factory, specs
from functools import partial
from lib import utils as utils


def bids_log(rule_name, **kwargs):
    """
    Args:
        rule_name (str): The name of the rule for the log.
        **kwargs: Optional parameters (e.g., subj_wildcards, hemi, label, etc.).

    Returns:
        str: Path to the log file.
    """

    log_params = {
        "root": "logs",
        "datatype": rule_name,
        "suffix": f"log.txt",
    }
    log_params.update(kwargs)

    return bids(**log_params)


def conda_env(env_name):
    """
    Returns the path to the Conda environment file if --use-conda is set, otherwise returns None.

    Args:
        env_name (str): The name of the environment.

    Returns:
        str or None: Path to the Conda YAML file or None.
    """
    import snakemake
    from packaging.version import parse

    # Get Snakemake version
    snakemake_version = parse(snakemake.__version__)

    if snakemake_version >= parse("8.0.0"):
        # Snakemake >= 8.0
        from snakemake.settings.types import DeploymentMethod

        if DeploymentMethod.CONDA in workflow.deployment_settings.deployment_method:
            return f"../envs/{env_name}.yaml"
    else:
        # Snakemake < 8.0
        if workflow.use_conda:
            return f"../envs/{env_name}.yaml"

    return None


def get_atlas_dir():
    return Path(utils.get_download_dir()) / "hippunfold-atlases"


def expand_hemi():
    if "hemi" in inputs[config["modality"]].zip_lists:
        # hemi is an input wildcard,
        #  so it will be already included when we expand
        return {}
    else:
        # hemi is not an input wildcard,
        # so we additionally expand using the config hemi
        return {"hemi": config["hemi"]}


def expand_label():
    return {"label": config["autotop_labels"]}


def get_single_bids_input(wildcards, component):
    """Generic input function for getting the first instance of a bids component
    from a dataset, and throw an error if more than one are found.

    Use this in a rule by using e.g.:
        in_img = partial(get_single_bids_input,component='T1w')
    """

    subj_inputs = inputs[component].filter(**wildcards).expand()
    if len(subj_inputs) > 1:
        raise ValueError(
            f"Expected 1 input for '{component}' {wildcards}, "
            f"but found {len(subj_inputs)}: {subj_inputs}."
            f"You can use the --filter-{component} option to filter your wildcards."
        )
    else:
        return subj_inputs[0]


def bids_atlas(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(
        Path(bids(root=root, tpl=template)) / bids(prefix=f"tpl-{template}", **entities)
    )


# take mean of all scans if >1, otherwise just copy the one scan
def get_avg_or_cp_scans_cmd(wildcards, input, output):
    if len(input) > 1:
        cmd = f"c3d {input} -mean -o {output}"
    else:
        cmd = f"cp {input} {output}"
    return cmd


def get_modality_suffix(modality):
    if modality[:4] == "hipp":
        return modality[4:]
    else:
        return modality


def get_inputs_spec_file(label, density):

    files = []
    files.extend(
        expand(
            bids(
                root=root,
                datatype="metric",
                den="{density}",
                suffix="{metric}.gii",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(label),
            density=density,
            allow_missing=True,
        )
    )
    files.extend(
        expand(
            bids(
                root=root,
                datatype="metric",
                den="{density}",
                desc="laplace",
                suffix="coords.shape.gii",
                hemi="{hemi}",
                dir="{dir}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            dir=["AP", "PD"],
            density=density,
            allow_missing=True,
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
            space=ref_spaces,
            density=density,
            allow_missing=True,
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
            surfname=["midthickness"],
            density=density,
            space=["unfold"],
            allow_missing=True,
        )
    )
    files.extend(
        expand(
            bids(
                root=root,
                datatype="cifti",
                den="{density}",
                suffix="{cifti}.nii",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(label),
            density=density,
            allow_missing=True,
        )
    )
    if label == "hipp":
        files.extend(
            expand(
                bids(
                    root=root,
                    datatype="cifti",
                    den="{density}",
                    atlas="{atlas}",
                    suffix="subfields.dlabel.nii",
                    label="{label}",
                    **inputs.subj_wildcards,
                ),
                atlas=config["atlas"],
                density=density,
                allow_missing=True,
            )
        )
        files.extend(
            expand(
                bids(
                    root=root,
                    datatype="metric",
                    den="{density}",
                    atlas="{atlas}",
                    suffix="subfields.label.gii",
                    hemi="{hemi}",
                    label="{label}",
                    **inputs.subj_wildcards,
                ),
                atlas=config["atlas"],
                density=density,
                allow_missing=True,
            )
        )

    return files


def get_final_spec():
    """include input files for spec too, so they don't get removed with temp()"""

    files = []

    files.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="surfaces.spec",
                **inputs.subj_wildcards,
            ),
            density=config["output_density"],
            allow_missing=True,
        )
    )
    # add spec inputs too
    for label in config["autotop_labels"]:
        for spec_input in get_inputs_spec_file(label, density=config["output_density"]):
            files.extend(
                inputs[config["modality"]].expand(
                    spec_input,
                    label=label,
                    **expand_hemi(),
                )
            )

    files.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="removeunused.touch",
                **inputs.subj_wildcards,
            )
        )
    )

    return files


def get_final_subfields():
    return inputs[config["modality"]].expand(
        bids(
            root=root,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="{space}",
            hemi="{hemi}",
            atlas="{atlas}",
            label="hipp",
            **inputs.subj_wildcards,
        ),
        hemi=config["hemi"],
        space=crop_ref_spaces,
        atlas=config["atlas"],
        allow_missing=True,
    )


def get_final_anat():
    anat = []
    if "T1w" in ref_spaces or "T2w" in ref_spaces:
        anat.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="anat",
                    desc="preproc",
                    suffix="{modality_suffix}.nii.gz".format(
                        modality_suffix=get_modality_suffix(config["modality"])
                    ),
                    space="{space}",
                    hemi="{hemi}",
                    **inputs.subj_wildcards,
                ),
                space=crop_ref_spaces,
                hemi=config["hemi"],
                allow_missing=True,
            )
        )
    return anat


def get_final_qc():
    qc = []

    if not template_modality == False:
        qc.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="qc",
                    suffix="regqc.png",
                    from_="{modality}",
                    to=config["template"],
                    **inputs.subj_wildcards,
                ),
                modality=template_modality,
                allow_missing=True,
            )
        )
    qc.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="qc",
                suffix="dseg.png",
                desc="subfields",
                space="{space}",
                hemi="{hemi}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            atlas=config["atlas"],
            space=crop_ref_spaces,
            allow_missing=True,
        )
    )
    qc.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="qc",
                suffix="midthickness.surf.png",
                den="{density}",
                desc="subfields",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            density=config["output_density"],
            label=config["autotop_labels"],
            space=ref_spaces,
            allow_missing=True,
        )
    )
    if len(config["hemi"]) == 2:
        qc.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="qc",
                    desc="subfields",
                    space="{space}",
                    atlas="{atlas}",
                    suffix="volumes.png",
                    **inputs.subj_wildcards,
                ),
                space=crop_ref_spaces,
                atlas=config["atlas"],
                allow_missing=True,
            )
        )
    if (config["modality"] == "T1w") or (config["modality"] == "T2w"):
        if not config["use_template_seg"]:
            qc.extend(
                inputs[config["modality"]].expand(
                    bids(
                        root=root,
                        datatype="qc",
                        desc="unetf3d",
                        suffix="dice.tsv",
                        hemi="{hemi}",
                        **inputs.subj_wildcards,
                    ),
                    hemi=config["hemi"],
                    allow_missing=True,
                )
            )
    return qc


def get_final_output():
    subj_output = []
    subj_output.extend(get_final_spec())
    subj_output.extend(get_final_subfields())
    subj_output.extend(get_final_anat())
    subj_output.extend(get_final_qc())
    return subj_output


if "corobl" in ref_spaces:

    ruleorder: laplace_beltrami > laynii_layers_equidist > laynii_layers_equivol


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


def get_create_atlas_output():
    files = str(
        Path(
            bids(
                root=get_atlas_dir(),
                tpl=config["new_atlas_name"],
            )
        )
        / "template_description.json"
    )
    return files


def get_input_for_shape_inject(wildcards):
    if config["modality"] == "dsegtissue":
        seg = bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    else:
        seg = bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    return seg


def get_cmd_warp_surface_2d_warp(wildcards, input, output):
    """Using this workaround for warping meshes with 2D warps, since surface-apply-warpfield was
    giving bounding box issues"""

    cmds = []
    cmds.append(
        f"wb_command -volume-to-surface-mapping {input.warp} {input.surf_gii} xywarp.shape.gii -trilinear"
    )
    cmds.append(
        f"wb_command -metric-math '0' zwarp.shape.gii -var DUMMY xywarp.shape.gii -column 1"
    )
    cmds.append(
        f"wb_command -metric-merge xyzwarp.shape.gii -metric xywarp.shape.gii  -metric zwarp.shape.gii"
    )
    cmds.append(
        f"wb_command -surface-coordinates-to-metric {input.surf_gii} coords.shape.gii"
    )
    cmds.append(
        f"wb_command -metric-math 'COORDS - WARP' warpedcoords.shape.gii -var COORDS coords.shape.gii -var WARP xyzwarp.shape.gii"
    )
    cmds.append(
        f"wb_command -surface-set-coordinates  {input.surf_gii} warpedcoords.shape.gii {output.surf_gii}"
    )
    return " && ".join(cmds)
