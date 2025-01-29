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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
            suffix="{modality,T1w|T2w}.nii.gz",
            space=config["template"],
            desc="affine"
        ),
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality,T1w|T2w}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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


# unfold ref nifti
rule create_unfold_ref:
    params:
        dims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["dims"]
        ),
        voxdims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["voxdims"]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["origin"]
        ),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.autotop]["orient"],
    output:
        nii=bids(
            root=work,
            space="unfold",
            label="{autotop}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}"

