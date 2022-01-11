
rule label_subfields_from_vol_coords_corobl:
    """ Label subfields using the volumetric coords and bigbrain labels"""
    input:
        subfields_mat=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "bigbrain",
            "BigBrain_ManualSubfieldsUnfolded.mat",
        ),
        nii_ap=bids(
            root=work,
            datatype="seg",
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
            datatype="seg",
            dir="PD",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    params:
        mat_name="subfields_avg",  #avg bigbrain over L/R hemis
    output:
        nii_label=bids(
            root=work,
            datatype="seg",
            desc="subfieldsnotissue",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    script:
        "../scripts/label_subfields_from_vol_coords.py"


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
            datatype="seg",
            desc="subfieldsnotissue",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    params:
        remap_dg="-threshold 8 8 6 0 -popas dg",
        remap_srlm="-threshold 2 2 7 0 -popas srlm",
        remap_cyst="-threshold 7 7 8 0 -popas cyst",
    output:
        combined=bids(
            root=work,
            datatype="seg",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
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
            datatype="seg",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        xfm=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=root,
            datatype="seg",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="{native_modality}",
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


rule resample_postproc_to_native:
    """Resample post-processed tissue seg to native"""
    input:
        nii=bids(
            root=work,
            datatype="seg",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}"
        ),
        xfm=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=work,
            datatype="seg",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="{native_modality}",
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
            datatype="seg",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}"
        ),
        xfm=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{native_modality}.nii.gz"
        ),
    output:
        nii=bids(
            root=work,
            datatype="seg",
            suffix="dseg.nii.gz",
            desc="unet",
            space="{native_modality}",
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


rule get_subfield_vols_subj:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            bids(
                root=root,
                **config["subj_wildcards"],
                datatype="seg",
                hemi="{hemi}",
                space="T2w",
                desc="subfields",
                suffix="dseg.nii.gz"
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
        lookup_tsv=os.path.join(
            workflow.basedir, "..", "resources", "desc-subfields_dseg.tsv"
        ),
    group:
        "subj"
    output:
        tsv=bids(
            root=root,
            datatype="seg",
            desc="subfields",
            suffix="volumes.tsv",
            **config["subj_wildcards"]
        ),
    script:
        "../scripts/gen_volume_tsv.py"


rule plot_subj_subfields:
    input:
        tsv=bids(
            root=root,
            datatype="seg",
            desc="subfields",
            suffix="volumes.tsv",
            **config["subj_wildcards"]
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                desc="subfields",
                suffix="volumes.png",
                **config["subj_wildcards"]
            ),
            caption="../report/subj_volume_plot.rst",
            category="Subfield Volumes",
        ),
    group:
        "subj"
    script:
        "../scripts/plot_subj_subfields.py"


def get_bg_img_for_subfield_qc(wildcards):

    if config["modality"][:3] == "seg":
        bg_modality = config["modality"][3:]
    else:
        bg_modality = config["modality"]

    return bids(
        root=root,
        datatype="seg",
        desc="preproc",
        suffix=f"{bg_modality}.nii.gz",
        space="crop{native_modality}",
        hemi="{hemi}",
        **config["subj_wildcards"],
    )


rule qc_subfield:
    input:
        img=get_bg_img_for_subfield_qc,
        seg=bids(
            root=root,
            datatype="seg",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="crop{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="dseg.png",
                desc="subfields",
                space="crop{native_modality}",
                hemi="{hemi}",
                **config["subj_wildcards"]
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    script:
        "../scripts/vis_qc_dseg.py"


rule qc_subfield_surf:
    input:
        surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            den="{density}",
            space="{native_modality}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="midthickness.surf.png",
                den="{density}",
                desc="subfields",
                space="crop{native_modality}",
                hemi="{hemi}",
                label="{autotop}",
                **config["subj_wildcards"]
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    script:
        "../scripts/vis_qc_surf.py"


rule concat_subj_vols_tsv:
    """Concatenate all subject tsv files into a single tsv"""
    input:
        tsv=lambda wildcards: expand(
            bids(
                root=root,
                datatype="seg",
                desc="subfields",
                suffix="volumes.tsv",
                **config["subj_wildcards"]
            ),
            subject=config["input_lists"][get_modality_key(config["modality"])][
                "subject"
            ],
            session=config["sessions"],
        ),
    group:
        "aggregate"
    output:
        tsv=bids(root=root, prefix="group", desc="subfields", suffix="volumes.tsv"),
    run:
        import pandas as pd

        pd.concat([pd.read_table(in_tsv) for in_tsv in input]).to_csv(
            output.tsv, sep="\t", index=False
        )


rule resample_subfields_to_unfold:
    """Resampling to unfold space"""
    input:
        nii=bids(
            root=work,
            datatype="seg",
            desc="subfields",
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        xfm=bids(
            root=work,
            datatype="seg",
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
            datatype="seg",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="unfold",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.xfm}  -t {input.xfm}"
