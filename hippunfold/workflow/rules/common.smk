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


def get_modality_key(modality):
    if modality[:3] == "seg":
        return "seg"
    else:
        return modality


def get_modality_suffix(modality):
    if modality[:3] == "seg":
        return modality[3:]
    elif modality[:4] == "hipp":
        return modality[4:]
    else:
        return modality


def get_final_spec():
    specs = inputs[get_modality_key(config["modality"])].expand(
        bids(
            root=root,
            datatype="surf",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            suffix="surfaces.spec",
            **inputs.subj_wildcards,
        ),
        density=config["output_density"],
        space=ref_spaces,
        hemi=config["hemi"],
        autotop=config["autotop_labels"],
        surfname=config["surf_types"]["hipp"],
        allow_missing=True,
    )
    specs.extent(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{autotop}",
                suffix="surfaces.spec",
                **inputs.subj_wildcards,
            ),
            density=config["output_density"],
            space=ref_spaces,
            hemi=config["hemi"],
            autotop=config["autotop_labels"],
            surfname=config["surf_types"]["dentate"],
            allow_missing=True,
        )
    )
    return specs


def get_final_subfields():
    return inputs[get_modality_key(config["modality"])].expand(
        bids(
            root=root,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="{space}",
            hemi="{hemi}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
        hemi=config["hemi"],
        space=crop_ref_spaces,
        atlas=config["atlas"],
        allow_missing=True,
    )


def get_final_coords():
    if "laplace" in config["laminar_coords_method"]:
        desc_io = "laplace"
    elif "equivolume" in config["laminar_coords_method"]:
        desc_io = "equivol"

    coords = []
    # compute all laplace coords by default (incl IO)
    coords.extend(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="coords",
                dir="{dir}",
                suffix="coords.nii.gz",
                desc="{desc}",
                space="{space}",
                hemi="{hemi}",
                label="{autotop}",
                **inputs.subj_wildcards,
            ),
            desc="laplace",
            dir=["AP", "PD", "IO"],
            autotop=config["autotop_labels"],
            hemi=config["hemi"],
            space=crop_ref_spaces,
            allow_missing=True,
        )
    )
    coords.extend(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="coords",
                dir="{dir}",
                suffix="coords.nii.gz",
                desc="{desc}",
                space="{space}",
                hemi="{hemi}",
                label="hipp",
                **inputs.subj_wildcards,
            ),
            desc=[desc_io],
            dir=["IO"],
            hemi=config["hemi"],
            space=crop_ref_spaces,
            allow_missing=True,
        )
    )
    return coords


def get_final_transforms():
    xfms = []

    xfms.extend(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                label="{autotop}",
                suffix="xfm.nii.gz",
                hemi="{hemi}",
                from_="{space}",
                to="unfold",
                mode="image",
            ),
            space=ref_spaces,
            autotop=config["autotop_labels"],
            hemi=config["hemi"],
            allow_missing=True,
        )
    )

    xfms.extend(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                label="{autotop}",
                suffix="xfm.nii.gz",
                hemi="{hemi}",
                from_="unfold",
                to="{space}",
                mode="image",
            ),
            space=ref_spaces,
            autotop=config["autotop_labels"],
            hemi=config["hemi"],
            allow_missing=True,
        )
    )

    xfms.extend(
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                label="{autotop}",
                suffix="refvol.nii.gz",
                space="unfold",
            ),
            autotop=config["autotop_labels"],
            allow_missing=True,
        )
    )

    return xfms


def get_final_anat():
    anat = []

    if "T1w" in ref_spaces or "T2w" in ref_spaces:
        anat.extend(
            inputs[get_modality_key(config["modality"])].expand(
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
            inputs[get_modality_key(config["modality"])].expand(
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
        inputs[get_modality_key(config["modality"])].expand(
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
        inputs[get_modality_key(config["modality"])].expand(
            bids(
                root=root,
                datatype="qc",
                suffix="midthickness.surf.png",
                den="{density}",
                desc="subfields",
                space="{space}",
                hemi="{hemi}",
                label="{autotop}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            autotop=config["autotop_labels"],
            density=config["output_density"],
            space=ref_spaces,
            allow_missing=True,
        )
    )
    if len(config["hemi"]) == 2:
        qc.extend(
            inputs[get_modality_key(config["modality"])].expand(
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
                inputs[get_modality_key(config["modality"])].expand(
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


def get_final_subj_output():
    subj_output = []
    subj_output.extend(get_final_spec())
    subj_output.extend(get_final_surf())
    subj_output.extend(get_final_subfields())
    subj_output.extend(get_final_coords())
    subj_output.extend(get_final_transforms())
    subj_output.extend(get_final_anat())
    subj_output.extend(get_final_qc())
    return subj_output


def get_final_output():
    if config["keep_work"]:
        subj_output = get_final_subj_output()
    else:
        subj_output = get_final_work_tar()

    final_output = []

    modality_suffix = get_modality_suffix(config["modality"])
    modality_key = get_modality_key(config["modality"])

    # use a zip list for subject/session:
    zip_list = inputs[modality_key].zip_lists
    if "session" in zip_list:
        zip_list = snakebids.filter_list(
            zip_list,
            {
                "session": inputs[get_modality_key(config["modality"])].zip_lists[
                    "session"
                ]
            },
        )

    final_output.extend(
        expand(
            expand(subj_output, zip, allow_missing=True, **zip_list),
            modality_suffix=modality_suffix,
        )
    )

    return final_output


if "corobl" in ref_spaces:

    rule copy_coords_to_results:
        input:
            os.path.join(work, "{pre}_space-corobl_{post}{suffix}.{ext}"),
        output:
            os.path.join(root, "{pre}_space-corobl_{post}{suffix,coords}.{ext}"),
        group:
            "subj"
        shell:
            "cp {input} {output}"

    rule copy_xfm_to_results:
        input:
            os.path.join(work, "{pre}_{fromto,from|to}-corobl_{post}{suffix}.{ext}"),
        output:
            os.path.join(root, "{pre}_{fromto,from|to}-corobl_{post}{suffix,xfm}.{ext}"),
        group:
            "subj"
        shell:
            "cp {input} {output}"

    rule copy_subfields_to_results:
        input:
            os.path.join(work, "{pre}_desc-subfields_{post}{suffix}.{ext}"),
        output:
            os.path.join(root, "{pre}_desc-subfields_{post}{suffix,dseg}.{ext}"),
        group:
            "subj"
        shell:
            "cp {input} {output}"


def get_final_work_tar():
    bids = bids_factory(specs.v0_0_0(subject_dir=False, session_dir=False))
    return bids(root=work, suffix="work.tar.gz", **inputs.subj_wildcards)


def get_work_dir(wildcards):
    folder_with_file = inputs[get_modality_key(config["modality"])].expand(
        bids(root=work, **inputs.subj_wildcards), **wildcards
    )
    folder_without_file = os.path.dirname(folder_with_file[0])
    return folder_without_file


def get_download_dir():
    if "HIPPUNFOLD_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["HIPPUNFOLD_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("hippunfold", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir


rule archive_work_after_final:
    input:
        get_final_subj_output(),
    params:
        work_dir=get_work_dir,
    output:
        get_final_work_tar(),
    group:
        "subj"
    shell:
        #exit code 0 or 1 is acceptable (2 is fatal)
        "tar -czf {output} {params.work_dir}; "
        "if [ $? -le 1 ]; then "
        "  rm -rf {params.work_dir}; "
        "else exit 1; "
        "fi"
