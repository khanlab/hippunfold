rule create_crop_ref:
    """Create ref space for hires crop in original modality space
    TODO:  expose the resampling factor and size as cmd line args"""
    input:
        seg=inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="subfields",
                space="{modality}",
                hemi="{hemi}",
                atlas="{atlas}",
                label="hipp",
                **inputs.subj_wildcards,
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
    params:
        resample=config["crop_res"],
        pad_to=config["crop_box"],
    output:
        ref=temp(
            bids(
                root=root,
                datatype="warps",
                suffix="cropref.nii.gz",
                space="{modality}",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -binarize -interpolation NearestNeighbor -trim 0vox -resample-mm {params.resample} -pad-to {params.pad_to} 0 -scale 0 -type uchar {output}"


rule resample_unet_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="unet",
                space="crop{modality}",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_postproc_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
            label=config["autotop_labels"][-1],
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="postproc",
                space="crop{modality}",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_subfields_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="crop{modality}",
            hemi="{hemi}",
            atlas="{atlas}",
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_coords_crop:
    input:
        nii=bids(
            root=root,
            datatype="coords",
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="coords",
            dir="{dir}",
            suffix="coords.nii.gz",
            desc="{desc}",
            space="crop{modality}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]"


rule resample_to_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{modality}.nii.gz",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="{modality}.nii.gz",
            space="crop{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/ants.yaml"
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
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T2w",
            to="{modality}",
            desc="rigid",
            type_="itk",
        )
    return xfm


rule resample_t2_to_crop:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="T2w.nii.gz",
            desc="preproc",
        ),
        ref=bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        xfm=get_xfm_t2_to_t1(),
    params:
        xfm_opt=lambda wildcards, input: (
            "" if len(input.xfm) == 0 else f"-t {input.xfm}"
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T2w.nii.gz",
            space="crop{modality}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output.nii} -r {input.ref} {params.xfm_opt}"
