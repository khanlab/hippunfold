

rule import_t1:
    """TODO: add motion-corrected averaging, like we do for T2w"""
    input:
        in_img=partial(get_single_bids_input, component="T1w"),
    output:
        temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T1w.nii.gz",
            )
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


if config["skip_preproc"]:

    rule import_preproc_t1:
        input:
            in_img=partial(get_single_bids_input, component="T1w"),
        output:
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T1w.nii.gz",
                desc="preproc",
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"

else:

    rule n4_t1:
        input:
            t1=bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T1w.nii.gz",
            ),
        output:
            t1=bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                desc="preproc",
                suffix="T1w.nii.gz",
            ),
        threads: 8
        container:
            config["singularity"]["autotop"]
        conda:
            conda_env("ants")
        group:
            "subj"
        shell:
            "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
            "N4BiasFieldCorrection -d 3 -i {input.t1} -o {output}"


# apply transform to get subject in corobl cropped space
rule warp_t1_to_corobl_crop:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz",
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["crop_ref"].format(
            **wildcards, modality="T1w"
        ),
    output:
        t1=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T1w.nii.gz",
                space="corobl",
                desc="preproc",
                hemi="{hemi,L|R}",
            )
        ),
    conda:
        conda_env("ants")
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.t1} -o {output.t1} -r {params.ref}  -t {input.xfm}"
