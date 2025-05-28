rule import_t2:
    input:
        inputs["T2w"].path,
    output:
        temp(
            bids(
                root=root,
                datatype="anat",
                suffix="T2w.nii.gz",
                **inputs["T2w"].wildcards,
            )
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule n4_t2:
    input:
        bids(
            root=root,
            datatype="anat",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
    output:
        temp(
            bids(
                root=root,
                datatype="anat",
                desc="n4",
                suffix="T2w.nii.gz",
                **inputs["T2w"].wildcards,
            )
        ),
    threads: 8
    conda:
        conda_env("ants")
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input} -o {output}"


def get_ref_n4_t2(wildcards):
    # get the first image
    t2_imgs = (
        inputs["T2w"]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                datatype="anat",
                desc="n4",
                suffix="T2w.nii.gz",
                **inputs["T2w"].wildcards,
            )
        )
    )

    return t2_imgs[0]


def get_floating_n4_t2(wildcards):
    t2_imgs = (
        inputs["T2w"]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                datatype="anat",
                suffix="T2w.nii.gz",
                desc="n4",
                **inputs["T2w"].wildcards,
            ),
        )
    )
    return t2_imgs[int(wildcards.idx)]


rule reg_t2_to_ref:
    input:
        ref=get_ref_n4_t2,
        flo=get_floating_n4_t2,
    output:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T2w{idx}",
            to="T2w0",
            desc="rigid",
            type_="ras",
        ),
        warped=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="T2w.nii.gz",
                desc="aligned",
                floating="{idx}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        conda_env("niftyreg")
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac"


rule ras_to_itk_reg_t2:
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T2w{idx}",
            to="T2w0",
            desc="rigid",
            type_="ras",
        ),
    output:
        xfm_itk=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w{idx}",
                to="T2w0",
                desc="rigid",
                type_="itk",
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d_affine_tool  {input.xfm_ras} -oitk {output.xfm_itk}"


def get_aligned_n4_t2(wildcards):
    # first get the number of floating t2s
    num_scans = len(inputs["T2w"].filter(**wildcards).expand())

    # then, return the path, expanding over range(1,num_scans) -i.e excludes 0 (ref image)
    t2_imgs = (
        inputs["T2w"]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                datatype="anat",
                desc="aligned",
                floating="{idx}",
                suffix="T2w.nii.gz",
                **inputs.subj_wildcards,
            ),
            idx=range(1, num_scans),
        )
    )
    return t2_imgs


if config["skip_preproc"]:

    # for preproc T2, we expect the user to use bids filters/wildcards to ensure only 1 T1w is matched per subject
    rule import_preproc_t2:
        input:
            lambda wildcards: inputs["T2w"].filter(**wildcards).expand()[0],
        output:
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T2w.nii.gz",
                desc="preproc",
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
                **inputs.subj_wildcards,
                suffix="T2w.nii.gz",
                desc="preproc",
            ),
        container:
            config["singularity"]["autotop"]
        conda:
            conda_env("c3d")
        group:
            "subj"
        shell:
            "{params.cmd}"


rule reg_t2_to_t1_part1:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="T2w.nii.gz",
            desc="preproc",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz",
        ),
    output:
        warped=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="T2w.nii.gz",
            desc="preproc",
            space="T1w",
        ),
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w",
                to="T1w",
                desc="rigid",
                type_="ras",
            )
        ),
    log:
        bids_log(
            "reg_t2_to_t1_part1",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("niftyreg")
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac &> {log}"


rule reg_t2_to_t1_part2:
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T2w",
            to="T1w",
            desc="rigid",
            type_="ras",
        ),
    output:
        xfm_itk=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w",
                to="T1w",
                desc="rigid",
                type_="itk",
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


def get_inputs_compose_t2_xfm_corobl(wildcards):
    if config["t1_reg_template"]:
        # xfm0: t2 to t1
        # xfm1: t1 to corobl
        t2_to_t1 = (
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w",
                to="T1w",
                desc="rigid",
                type_="itk",
            ),
        )
        t1_to_cor = (
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
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
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
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
        t2_to_cor=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w",
                to="corobl",
                desc="affine",
                type_="itk",
            )
        ),
    log:
        bids_log(
            "compose_t2_xfm_corobol",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "{params.cmd} > {log}"


# if already have t2w in T1w space, then we don't need to use composed xfm:
def get_xfm_to_corobl():
    if config["skip_coreg"]:
        xfm = bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk",
        )
    else:
        xfm = (
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="T2w",
                to="corobl",
                desc="affine",
                type_="itk",
            ),
        )
    return xfm


# apply transform to get subject in corobl cropped space -- note that this could be marginally improved if we xfm each T2 to the cropped space in single resample (instead of avgT2w first)
rule warp_t2_to_corobl_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="T2w.nii.gz",
            desc="preproc",
        ),
        xfm=get_xfm_to_corobl(),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["crop_ref"].format(**wildcards),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="T2w.nii.gz",
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
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output.nii} -r {params.ref}  -t {input.xfm}"
