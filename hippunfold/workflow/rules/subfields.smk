rule import_dseg_subfields:
    """read in volumetric subfield dseg, used for group_create_atlas only"""
    input:
        vol_dseg=partial(get_single_bids_input, component="dsegsubfields"),
        label_list=Path(workflow.basedir)
        / "../resources/label_lut/labellist_withdg.txt",
    output:
        label_dseg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                desc="subfields",
                suffix="dseg.nii.gz",
                space="corobl",
                hemi="{hemi,L|R}",
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -volume-label-import {input.vol_dseg} {input.label_list} {output.label_dseg}"


rule subfields_to_label_gifti:
    input:
        vol=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
        ),
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    output:
        label_gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="subfields.label.gii",
                den="native",
                hemi="{hemi}",
                label="{label,hipp}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -volume-label-to-surface-mapping {input.vol} {input.surf_gii} {output.label_gii}"


rule native_label_gii_to_unfold_nii:
    """converts subfields label .gii files to .nii for applying with 2D ANTS transform (when creating a template)"""
    input:
        label_gii=bids(
            root=root,
            datatype="metric",
            suffix="subfields.label.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_nii=bids(
            root=root,
            datatype="warps",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        interp="-nearest-vertex 10",
    output:
        label_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="subfieldsfromnative.nii.gz",
                space="unfold",
                hemi="{hemi}",
                label="{label,hipp}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {input.label_gii} {input.midthickness_surf} {input.ref_nii} {output.label_nii} "
        " {params.interp}"


rule unfoldreg_label_gii_to_unfold_nii:
    """converts subfields label .gii files from atlas to subject unfold space, for creating template with new_atlas_subfields_from==unfoldreg"""
    input:
        label_gii=bids(
            root=root,
            datatype="metric",
            suffix="subfields.label.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            atlas=config["atlas"],
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="unfoldreg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfoldreg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="unfoldreg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_nii=bids(
            root=root,
            datatype="warps",
            suffix="refvol.nii.gz",
            space="unfold",
            desc="slice",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        interp="-nearest-vertex 10",
    output:
        label_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="subfieldsfromunfoldreg.nii.gz",
                space="unfold",
                hemi="{hemi}",
                label="{label,hipp}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {input.label_gii} {input.midthickness_surf} {input.ref_nii} {output.label_nii} "
        " {params.interp}"


rule label_subfields_from_vol_coords_corobl:
    input:
        ref_nii=get_labels_for_laplace,
        midthickness_surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        inner_surf=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer_surf=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        label_gii=bids(
            root=root,
            datatype="metric",
            suffix="subfields.label.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    output:
        nii_label=temp(
            bids(
                root=root,
                datatype="anat",
                desc="subfieldsnotissue",
                suffix="dseg.nii.gz",
                space="corobl",
                hemi="{hemi}",
                atlas="{atlas}",
                label="{label,hipp}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    log:
        bids_log(
            "label_subfields_from_vol_coords_corobl",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
        ),
    shell:
        "wb_command -label-to-volume-mapping {input.label_gii} {input.midthickness_surf} {input.ref_nii} {output.nii_label} &>> {log}"
        " -ribbon-constrained {input.inner_surf} {input.outer_surf} &>> {log}"


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
            root=root,
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
        combined=temp(
            bids(
                root=root,
                datatype="anat",
                desc="subfields",
                suffix="dseg.nii.gz",
                space="corobl",
                hemi="{hemi}",
                label="{label,hipp}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input.tissue} -dup -dup {params.remap} {input.subfields} -push dg -max -push srlm -max -push cyst -max -type uchar -o {output}"


rule resample_subfields_to_orig:
    """Resampling to original modality space"""
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
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{modality}.nii.gz",
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="{modality,T1w|T2w}",
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


rule resample_postproc_to_orig:
    """Resample post-processed tissue seg to original modality"""
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
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{modality}.nii.gz",
        ),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="postproc",
                space="{modality,T2w|T2w}",
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


rule resample_unet_to_orig:
    """Resample unet tissue seg to original modality"""
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
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{modality}.nii.gz",
        ),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="unet",
                space="{modality,T1w|T2w}",
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


rule resample_subfields_to_unfold:
    """Resampling to unfold space"""
    input:
        nii=bids(
            root=root,
            datatype="anat",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
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
        "../envs/ants.yaml"
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.xfm}  -t {input.xfm}"
