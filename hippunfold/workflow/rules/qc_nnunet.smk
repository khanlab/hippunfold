"""Quality control rules for nnUNet segmentation"""
from pathlib import Path


def get_f3d_ref(wildcards, input):
    """Get reference image for F3D registration based on modality"""
    if config["modality"] == "T2w":
        nii = Path(input.template_dir) / config["template_files"][config["template"]][
            "crop_ref"
        ].format(**wildcards)
    elif config["modality"] == "T1w":
        nii = Path(input.template_dir) / config["template_files"][config["template"]][
            "crop_refT1w"
        ].format(**wildcards)
    else:
        raise ValueError("modality not supported for nnunet!")
    return nii


rule qc_nnunet_f3d:
    """Register nnUNet segmentation to template space for quality control"""
    input:
        img=(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="{modality}.nii.gz".format(modality=config["modality"]),
                space="corobl",
                desc="preproc",
                hemi="{hemi}",
            ),
        ),
        seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=get_f3d_ref,
    output:
        cpp=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="cpp.nii.gz",
                desc="f3d",
                space="corobl",
                hemi="{hemi}",
            )
        ),
        res=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="{modality}.nii.gz".format(modality=config["modality"]),
                desc="f3d",
                space="template",
                hemi="{hemi}",
            )
        ),
        res_mask=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="mask.nii.gz",
                desc="f3d",
                space="template",
                hemi="{hemi}",
            )
        ),
    conda:
        "../envs/niftyreg.yaml"
    log:
        bids_log(
            "qc_nnunet_f3d",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    group:
        "subj"
    shell:
        "reg_f3d -flo {input.img} -ref {params.ref} -res {output.res} -cpp {output.cpp} &> {log} && "
        "reg_resample -flo {input.seg} -cpp {output.cpp} -ref {params.ref} -res {output.res_mask} -inter 0 &> {log}"


rule qc_nnunet_dice:
    """Calculate Dice coefficient between nnUNet segmentation and template"""
    input:
        res_mask=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="mask.nii.gz",
            desc="f3d",
            space="template",
            hemi="{hemi}",
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        hipp_lbls=[1, 2, 7, 8],
        ref=lambda wildcards, input: (
            Path(input.template_dir)
            / config["template_files"][config["template"]]["Mask_crop"].format(
                **wildcards
            )
        ),
    output:
        dice=report(
            bids(
                root=root,
                datatype="qc",
                suffix="dice.tsv",
                desc="unetf3d",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            ),
            caption="../report/nnunet_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    conda:
        "../envs/pyunfold.yaml"
    script:
        "../scripts/dice.py"
