rule import_dseg_subfields:
    """read in volumetric subfield dseg, used for participant_create_template only"""
    input:
        vol_dseg=partial(get_single_bids_input, component="dsegsubfields"),
        label_list=Path(workflow.basedir) / "../resources/atlas-v2/labellist_withdg.txt",
    output:
        label_dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "wb_command -volume-label-import {input.vol_dseg} {input.label_list} {output.label_dseg}"


rule subfields_to_label_gifti:
    input:
        vol=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
        ),
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("workbench")
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -volume-label-to-surface-mapping {input.vol} {input.surf_gii} {output.label_gii}"


rule label_subfields_from_vol_coords_corobl:
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_nii=get_labels_for_laplace,
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    output:
        nii_label=bids(
            root=work,
            datatype="anat",
            desc="subfieldsnotissue",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {input.label_gii} {input.midthickness_surf} {input.ref_nii} {output.nii_label} "
        " -ribbon-constrained {input.inner_surf} {input.outer_surf}"


def get_tissue_atlas_remapping(wildcards):
    mapping = config["tissue_atlas_mapping"]

    remap = []

    for label in mapping["tissue"].keys():
        in_label = mapping["tissue"][label]
        out_label = mapping.get(wildcards.atlas, mapping.get("default"))[label]

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
            label="{label}",
            **inputs.subj_wildcards,
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
            label="{label,hipp}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("c3d")
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
            label="{label}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{native_modality}.nii.gz",
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
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("ants")
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
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
            label=config["autotop_labels"][-1],
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{native_modality}.nii.gz",
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="{native_modality,T2w|T2w}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("ants")
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
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{native_modality}.nii.gz",
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="unet",
            space="{native_modality,T1w|T2w}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("ants")
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
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="image",
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
            **inputs.subj_wildcards,
        ),

    conda:
        conda_env("ants")
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.xfm}  -t {input.xfm}"
