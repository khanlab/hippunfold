rule create_native_crop_ref:
    """Create ref space for hires crop in native space
    TODO:  expose the resampling factor and size as cmd line args"""
    input:
        seg=expand(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="subfields",
                space="{native_modality}",
                hemi="{hemi}",
                atlas="{atlas}",
                **config["subj_wildcards"]
            ),
            atlas=config["atlas"][0],
            allow_missing=True,
        ),
    params:
        resample=config["crop_native_res"],
        pad_to=config["crop_native_box"],
    output:
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -binarize -interpolation NearestNeighbor -trim 0vox -resample-mm {params.resample} -pad-to {params.pad_to} 0 -scale 0 -type uchar {output}"


rule resample_unet_native_crop:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}"
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="unet",
            space="crop{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_postproc_native_crop:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}"
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="crop{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_subfields_native_crop:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="crop{native_modality}",
            hemi="{hemi}",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_coords_native_crop:
    input:
        nii=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=root,
            datatype="coords",
            dir="{dir}",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="crop{native_modality}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_native_to_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="{native_modality}.nii.gz"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="{native_modality}.nii.gz",
            space="crop{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output.nii} -r {input.ref} "


def get_xfm_t2_to_t1():
    if config["skip_coreg"]:
        xfm = []
    else:
        xfm = bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="T2w",
            to="{native_modality}",
            desc="rigid",
            type_="itk"
        )
    return xfm


rule resample_t2_to_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="T2w.nii.gz",
            desc="preproc"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        xfm=get_xfm_t2_to_t1(),
    params:
        xfm_opt=lambda wildcards, input: ""
        if len(input.xfm) == 0
        else f"-t {input.xfm}",
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T2w.nii.gz",
            space="crop{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output.nii} -r {input.ref} {params.xfm_opt}"
