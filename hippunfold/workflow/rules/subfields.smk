
rule resample_unfoldreg_subfields:
    input:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="{atlas}",
            to="subject",
            space="unfold",
            type_="itk",
            hemi="{hemi}"
        ),
    output:
        label_nii=bids(
            root=root,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        "c3d {input.label_nii} -slice z 0:15 -oo tmp0%d.nii.gz && "
        "for fn in $(ls tmp*.nii.gz); do antsApplyTransforms -d 2 -i $fn -r $fn -o $fn -n MultiLabel -t {input.warp}; done && "
        "c3d tmp*.nii.gz -tile z -o recombined.nii.gz && "
        "c3d {input.label_nii} recombined.nii.gz -copy-transform -o {output}"


def skip_unfoldreg_option_subfields(wildcards):
    if config["no_unfolded_reg"]:
        label_nii = bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        )
    else:
        label_nii = bids(
            root=root,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            atlas="{atlas}",
            **config["subj_wildcards"]
        )
    return label_nii


rule label_subfields_from_vol_coords_corobl:
    """ Label subfields using the volumetric coords and atlas subfield labels"""
    input:
        label_nii=skip_unfoldreg_option_subfields,
        nii_ap=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        nii_pd=bids(
            root=work,
            datatype="coords",
            dir="PD",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        nii_io=get_laminar_coords,
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: config["laplace_labels"]["AP"]["gm"],
    output:
        nii_label=bids(
            root=work,
            datatype="anat",
            desc="subfieldsnotissue",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/label_subfields_from_vol_coords.py"


def get_tissue_atlas_remapping(wildcards):
    mapping = config["tissue_atlas_mapping"]

    remap = []

    for label in mapping["tissue"].keys():
        in_label = mapping["tissue"][label]
        out_label = mapping[wildcards.atlas][label]

        remap.append(f"-threshold {in_label} {in_label} {out_label} 0 -popas {label}")

    return " ".join(remap)


rule combine_tissue_subfield_labels_corobl:
    """Combine subfield labels with the DG, SRLM and Cyst tissue labels

    add srlm, cyst, dg from postproc labels to subfields

    first remap tissue labels to get three sep labels
    then, we just need to add those in, using max(old,new) to override old with new in case of conflict
    """
    input:
        tissue=get_labels_for_laplace,
        subfields=bids(
            root=work,
            datatype="anat",
            desc="subfieldsnotissue",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    params:
        remap=get_tissue_atlas_remapping,
    output:
        combined=bids(
            root=work,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.tissue} -dup -dup {params.remap} {input.subfields} -push dg -max -push srlm -max -push cyst -max -type uchar -o {output}"


rule resample_subfields_to_native:
    """Resampling to native space"""
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
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="{native_modality,T1w|T2w}",
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


rule resample_postproc_to_native:
    """Resample post-processed tissue seg to native"""
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
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="{native_modality,T2w|T2w}",
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


rule resample_unet_to_native:
    """Resample unet tissue seg to native"""
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
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            desc="preproc",
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="unet",
            space="{native_modality,T1w|T2w}",
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


rule resample_subfields_to_unfold:
    """Resampling to unfold space"""
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
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="image"
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="unfold",
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
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.xfm}  -t {input.xfm}"
