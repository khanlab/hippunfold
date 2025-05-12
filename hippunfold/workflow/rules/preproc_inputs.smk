rule import_any_modality:
    input:
        inputs[config["modality"]].path,
    output:
        bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        )
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule lamareg_to_template:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["T1w"].format(**wildcards),
    output:
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.txt",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invaffine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.txt",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
    group:
        "subj" 
    shell:
        "lamar generate-warpfield --fixed {params.ref} --moving {input.img} --affine-file {output.affine} --rev-affine-file {output.invaffine} --warpfield {output.warp} --ref-warp-file {output.invwarp}"


rule apply_transforms:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.txt",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["Mask_crop"].format(**wildcards),
    output:
        img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    group:
        "subj"
    shell:
        "lamar apply-warp --affine-file {input.affine} --moving {input.img} --reference {params.ref} --output-file {output.img}"


# TODO: refine registrations using the raw images and the above initializations


rule template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            type_="itk",
        ),
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="{modality,T1w|T2w}",
                to="corobl",
                type_="ras",
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -o {output}"


rule superres_inputs:
    input:
        img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    output:
        img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix="preproc.nii.gz",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input.img} -add -o {output.img}"
