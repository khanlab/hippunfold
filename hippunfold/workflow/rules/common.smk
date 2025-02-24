from appdirs import AppDirs
from snakebids.paths import bids_factory, specs

from functools import partial


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


def get_final_spec():
    specs = []
    specs.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                space="{space}",
                suffix="surfaces.spec",
                **inputs.subj_wildcards,
            ),
            density=config["output_density"],
            space=ref_spaces,
            allow_missing=True,
        )
    )
    specs.extend(
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                space="{space}",
                suffix="surfaces.spec",
                **inputs.subj_wildcards,
            ),
            space="corobl",
            allow_missing=True,
        )
    )

    return specs


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
                    from_="{native_modality}",
                    to=config["template"],
                    **inputs.subj_wildcards,
                ),
                native_modality=template_modality,
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
            label=config["autotop_labels"],
            density=config["output_density"],
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

    ruleorder: laplace_beltrami > laynii_layers_equidist > laynii_layers_equivol > copy_coords_to_results

    rule copy_coords_to_results:
        input:
            os.path.join(work, "{pre}_space-corobl_{post}{suffix}.{ext}"),
        output:
            os.path.join(root, "{pre,[^/].+}_space-corobl_{post}{suffix,coords}.{ext}"),
        group:
            "subj"
        shell:
            "cp {input} {output}"

    rule copy_xfm_to_results:
        input:
            os.path.join(work, "{pre}_{fromto}-corobl_{post}{suffix}.{ext}"),
        output:
            os.path.join(
                root, "{pre,[^/].+}_{fromto,from|to}-corobl_{post}{suffix,xfm}.{ext}"
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"

    rule copy_subfields_to_results:
        input:
            os.path.join(work, "{pre}_desc-subfields_{post}{suffix}.{ext}"),
        output:
            os.path.join(root, "{pre,[^/].+}_desc-subfields_{post}{suffix,dseg}.{ext}"),
        group:
            "subj"
        shell:
            "cp {input} {output}"


def get_download_dir():
    if "HIPPUNFOLD_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["HIPPUNFOLD_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("hippunfold", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir


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


def get_create_template_output():

    files = []
    for label in config["autotop_labels"]:
        files.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="surf",
                    suffix="{metric}.gii",
                    space="{space}",
                    hemi="{hemi_}",
                    label=label,
                    **inputs.subj_wildcards,
                ),
                metric=get_gifti_metric_types(label),
                space="corobl",
                hemi_=config["hemi"],
            )
        )
        files.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="surf",
                    suffix="{surftype}.surf.gii",
                    space="{space}",
                    hemi="{hemi_}",
                    label=label,
                    **inputs.subj_wildcards,
                ),
                metric=get_gifti_metric_types(label),
                space=["corobl", "unfold"],
                surftype=["inner", "outer", "midthickness"],
                hemi_=config["hemi"],
            )
        )
        files.extend(
            inputs[config["modality"]].expand(
                bids(
                    root=root,
                    datatype="surf",
                    suffix="subfields.label.gii",
                    space="corobl",
                    hemi="{hemi_}",
                    label="hipp",
                    **inputs.subj_wildcards,
                ),
                space="corobl",
                hemi_=config["hemi"],
            )
        )
        files.extend(
            expand(
            "template/pairs_{hemi}_{label}.csv",
            hemi=config["hemi"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
        )),
    return files


def get_input_for_shape_inject(wildcards):
    if config["modality"] == "dsegtissue":
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    return seg
