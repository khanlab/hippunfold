rule import_seg:
    input:
        lambda wildcards: expand(
            config["input_path"]["seg"],
            zip,
            **snakebids.filter_list(config["input_zip_lists"]["seg"], wildcards)
        )[0],
    output:
        bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["seg"],
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


# now, transform to coronal oblique, xfm dependent on the space wildcard
rule warp_seg_to_corobl_crop:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["input_wildcards"]["seg"],
            suffix="dseg.nii.gz"
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{space}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["crop_ref"].format(**wildcards),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi,L|R}",
            from_="{space}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {params.ref}  -t {input.xfm}"


rule lr_flip_seg:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            from_="{space}"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi,L}flip",
            from_="{space}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -flip x -o  {output}"
