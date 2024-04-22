rule import_t2:
    input:
        config["input_path"]["T2w"],
    output:
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["T2w"],
            suffix="T2w.nii.gz"
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule n4_t2:
    input:
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["T2w"],
            suffix="T2w.nii.gz"
        ),
    output:
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["T2w"],
            suffix="T2w.nii.gz",
            desc="n4"
        ),
    threads: 8
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input} -o {output}"


def get_ref_n4_t2(wildcards):
    # get the first image
    t2_imgs = expand(
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["T2w"],
            suffix="T2w.nii.gz",
            desc="n4"
        ),
        zip,
        **snakebids.filter_list(config["input_zip_lists"]["T2w"], wildcards)
    )
    return t2_imgs[0]


def get_floating_n4_t2(wildcards):
    t2_imgs = expand(
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["T2w"],
            suffix="T2w.nii.gz",
            desc="n4"
        ),
        zip,
        **snakebids.filter_list(config["input_zip_lists"]["T2w"], wildcards)
    )
    return t2_imgs[int(wildcards.idx)]


rule reg_t2_to_ref:
    input:
        ref=get_ref_n4_t2,
        flo=get_floating_n4_t2,
    output:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w{idx}",
            to="T2w0",
            desc="rigid",
            type_="ras"
        ),
        xfm_itk=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w{idx}",
            to="T2w0",
            desc="rigid",
            type_="itk"
        ),
        warped=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="aligned",
            floating="{idx}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac && "
        "c3d_affine_tool  {output.xfm_ras} -oitk {output.xfm_itk}"


def get_aligned_n4_t2(wildcards):
    # first get the number of floating t2s
    filtered = snakebids.filter_list(config["input_zip_lists"]["T2w"], wildcards)
    num_scans = len(filtered["subject"])

    # then, return the path, expanding over range(1,num_scans) -i.e excludes 0 (ref image)
    t2_imgs = expand(
        bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="aligned",
            floating="{idx}"
        ),
        **wildcards,
        idx=range(1, num_scans)
    )
    return t2_imgs


if config["skip_preproc"]:

    # grabs the first t2w only
    rule import_preproc_t2:
        input:
            lambda wildcards: expand(
                config["input_path"]["T2w"],
                zip,
                **snakebids.filter_list(config["input_zip_lists"]["T2w"], wildcards)
            )[0],
        output:
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="T2w.nii.gz",
                desc="preproc"
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"

else:

    # average aligned n4 images to get preproc T2w (or if single scan, just copy it)
    rule avg_aligned_or_cp_t2:
        input:
            ref=get_ref_n4_t2,
            flo=get_aligned_n4_t2,
        params:
            cmd=get_avg_or_cp_scans_cmd,
        output:
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="T2w.nii.gz",
                desc="preproc"
            ),
        container:
            config["singularity"]["autotop"]
        group:
            "subj"
        shell:
            "{params.cmd}"


rule reg_t2_to_t1:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="preproc"
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
    output:
        warped=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="preproc",
            space="T1w"
        ),
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w",
            to="T1w",
            desc="rigid",
            type_="ras"
        ),
        xfm_itk=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w",
            to="T1w",
            desc="rigid",
            type_="itk"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="reg.txt",
            from_="T2w",
            to="T1w",
            desc="rigid",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac &> {log} && "
        "c3d_affine_tool  {output.xfm_ras} -oitk {output.xfm_itk}"


def get_inputs_compose_t2_xfm_corobl(wildcards):
    if config["t1_reg_template"]:
        # xfm0: t2 to t1
        # xfm1: t1 to corobl
        t2_to_t1 = (
            bids(
                root=work,
                datatype="warps",
                **config["subj_wildcards"],
                suffix="xfm.txt",
                from_="T2w",
                to="T1w",
                desc="rigid",
                type_="itk",
            ),
        )
        t1_to_cor = (
            bids(
                root=work,
                datatype="warps",
                **config["subj_wildcards"],
                suffix="xfm.txt",
                from_="T1w",
                to="corobl",
                desc="affine",
                type_="itk",
            ),
        )
        return {"t2_to_t1": t2_to_t1, "t1_to_cor": t1_to_cor}

    else:
        # xfm0: t2 to template
        t2_to_std = (
            bids(
                root=work,
                datatype="warps",
                **config["subj_wildcards"],
                suffix="xfm.txt",
                from_="T2w",
                to=config["template"],
                desc="affine",
                type_="itk",
            ),
        )

        # xfm1: template to corobl
        template_dir = Path(download_dir) / "template" / config["template"]
        return {"t2_to_std": t2_to_std, "template_dir": template_dir}


def get_cmd_compose_t2_xfm_corobl(wildcards, input, output):
    if config["t1_reg_template"]:
        # xfm0: t2 to t1
        xfm0 = input.t2_to_t1
        # xfm1: t1 to corobl
        xfm1 = input.t1_to_cor
    else:
        # xfm0: t2 to template
        xfm0 = input.t2_to_std
        # xfm1: template to corobl
        xfm1 = Path(input.template_dir) / config["template_files"][config["template"]][
            "xfm_corobl"
        ].format(**wildcards)

    return f"c3d_affine_tool -itk {xfm0} -itk {xfm1} -mult -oitk {output}"


# now have t2 to t1 xfm, compose this with t1 to corobl xfm
rule compose_t2_xfm_corobl:
    input:
        unpack(get_inputs_compose_t2_xfm_corobl),
    params:
        cmd=get_cmd_compose_t2_xfm_corobl,
    output:
        t2_to_cor=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="composecorobl.txt",
            from_="T2w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd} > {log}"


# if already have t2w in T1w space, then we don't need to use composed xfm:
def get_xfm_to_corobl():
    if config["skip_coreg"]:
        xfm = bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        )
    else:
        xfm = (
            bids(
                root=work,
                datatype="warps",
                **config["subj_wildcards"],
                suffix="xfm.txt",
                from_="T2w",
                to="corobl",
                desc="affine",
                type_="itk"
            ),
        )
    return xfm


# apply transform to get subject in corobl cropped space -- note that this could be marginally improved if we xfm each T2 to the cropped space in single resample (instead of avgT2w first)
rule warp_t2_to_corobl_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="preproc"
        ),
        xfm=get_xfm_to_corobl(),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["crop_ref"].format(**wildcards),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
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
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output.nii} -r {params.ref}  -t {input.xfm}"


rule lr_flip_t2:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            space="corobl",
            desc="{desc}",
            hemi="{hemi}"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
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
