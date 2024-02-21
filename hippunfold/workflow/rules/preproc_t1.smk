ruleorder: compose_template_xfm_corobl > convert_template_xfm_ras2itk


rule import_t1:
    """Note: this rule only grabs the first T1w
    TODO: add motion-corrected averaging, like we do for T2w"""
    input:
        lambda wildcards: expand(
            config["input_path"]["T1w"],
            zip,
            **snakebids.filter_list(config["input_zip_lists"]["T1w"], wildcards)
        )[0],
    output:
        bids(
            root=work, datatype="anat", **config["subj_wildcards"], suffix="T1w.nii.gz"
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


if config["skip_preproc"]:

    rule import_preproc_t1:
        input:
            lambda wildcards: expand(
                config["input_path"]["T1w"],
                zip,
                **snakebids.filter_list(config["input_zip_lists"]["T1w"], wildcards)
            )[0],
        output:
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="T1w.nii.gz",
                desc="preproc"
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"

else:

    rule n4_t1:
        input:
            t1=bids(
                root=work,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="T1w.nii.gz"
            ),
        output:
            t1=bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                desc="preproc",
                suffix="T1w.nii.gz"
            ),
        threads: 8
        container:
            config["singularity"]["autotop"]
        group:
            "subj"
        shell:
            "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
            "N4BiasFieldCorrection -d 3 -i {input.t1} -o {output}"


def reg_to_template_cmd(wildcards, input, output):
    ref = str(
        Path(input.template_dir)
        / config["template_files"][config["template"]][wildcards.modality].format(
            **wildcards
        ),
    )
    if config["no_reg_template"]:
        cmd = f"reg_resample -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {input.xfm_identity}; cp {input.xfm_identity} {output.xfm_ras}"
    elif config["rigid_reg_template"]:
        cmd = f"reg_aladin -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {output.xfm_ras} -rigOnly"
    else:
        cmd = f"reg_aladin -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {output.xfm_ras}"
    return cmd


rule reg_to_template:
    """ generic for T1w or T2w right now """
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="{modality}.nii.gz"
        ),
        xfm_identity=os.path.join(workflow.basedir, "..", config["xfm_identity"]),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        cmd=reg_to_template_cmd,
    output:
        warped_subj=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{modality,T1w|T2w}.nii.gz",
            space=config["template"],
            desc="affine"
        ),
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{modality,T1w|T2w}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="reg.txt",
            from_="{modality,T1w|T2w}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


rule convert_template_xfm_ras2itk:
    input:
        bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{reg_suffix}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    output:
        bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{reg_suffix}",
            to=config["template"],
            desc="affine",
            type_="itk"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input}  -oitk {output}"


# now have subject -> template transform, can compose that with template -> corobl to get subject -> corobl
rule compose_template_xfm_corobl:
    input:
        sub_to_std=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to=config["template"],
            desc="affine",
            type_="itk"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        std_to_cor=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
    output:
        sub_to_cor=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input.sub_to_std} -itk {params.std_to_cor} -mult -oitk {output}"


rule invert_template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    output:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affineInverse",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -inv -o {output}"


rule template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    output:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality,T1w|T2w}",
            to="corobl",
            desc="affine",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -o {output}"


# apply transform to get subject in corobl cropped space
rule warp_t1_to_corobl_crop:
    input:
        t1=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["crop_ref"].format(
            **wildcards, modality="T1w"
        ),
    output:
        t1=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T1w.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi,L|R}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.t1} -o {output.t1} -r {params.ref}  -t {input.xfm}"


rule lr_flip_t1:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T1w.nii.gz",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T1w.nii.gz",
            space="corobl",
            desc="{desc}",
            hemi="{hemi,L}flip"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -flip x -o  {output}"
