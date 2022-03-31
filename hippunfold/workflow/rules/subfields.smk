

rule label_subfields_from_unfolded_atlas:
    input:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="hipp",
            **config["subj_wildcards"]
        ),
        warpitk_unfold2native=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            label="dentate",
            suffix="xfm.nii.gz",
            hemi="{hemi,Lflip|R}",
            from_="unfold",
            to="corobl",
            mode="image"
        ),
        refvol=get_labels_for_laplace,
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
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.label_nii} -r {input.refvol} -t {input.warpitk_unfold2native} -o {output.nii_label}"


rule combine_tissue_subfield_labels_corobl:
    """Combine subfield labels with the DG, SRLM and Cyst tissue labels

    add srlm, cyst, dg from postproc labels to subfields
    input dg label 8, output 6
    input srlm label 2, output 7
    input cyst label 7, output 8

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
        remap_dg="-threshold 8 8 6 0 -popas dg",
        remap_srlm="-threshold 2 2 7 0 -popas srlm",
        remap_cyst="-threshold 7 7 8 0 -popas cyst",
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
        "c3d {input.tissue} -dup {params.remap_dg} -dup {params.remap_srlm} {params.remap_cyst} {input.subfields} -push dg -max -push srlm -max -push cyst -max -type uchar -o {output}"


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
