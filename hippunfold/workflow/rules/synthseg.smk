"""SynthSeg v1 inference rules"""
import re
from pathlib import Path


def get_synthseg_input(wildcards):
    """Get the appropriate input image for SynthSeg inference"""
    
    # SynthSeg typically works with T1w or T2w images
    if config["modality"] == "T2w":
        return bids(
            root=root,
            datatype="anat",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="T2w.nii.gz",
            **inputs.subj_wildcards,
        ).format(**wildcards)
    elif config["modality"] == "T1w":
        return bids(
            root=root,
            datatype="anat",
            space="corobl",
            desc="preproc", 
            hemi="{hemi}",
            suffix="T1w.nii.gz",
            **inputs.subj_wildcards,
        ).format(**wildcards)
    elif config["modality"] == "hippb500":
        return bids(
            root=root,
            datatype="dwi",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="b500.nii.gz",
            **inputs.subj_wildcards,
        ).format(**wildcards)
    elif config["modality"] == "dsegtissue":
        return bids(
            root=root,
            datatype="anat",
            space="corobl",
            desc="cropped",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
            **inputs.subj_wildcards,
        ).format(**wildcards)
    else:
        raise ValueError(f"Unsupported modality for SynthSeg: {config['modality']}")


def get_cmd_copy_inputs_synthseg(wildcards, input):
    """Copy input images with SynthSeg naming convention"""
    in_img = input.in_img
    input_nii = "input.nii.gz"
    
    cmd = f"cp {in_img} {input_nii}"
    
    return cmd


rule run_inference_synthseg:
    """SynthSeg inference with checkpoint extraction in shadow directory"""
    input:
        in_img=get_synthseg_input,
        model_tar=get_model_tar(),
    params:
        model_dir="tempmodel",
        checkpoint_path="tempmodel/synthseg/last.ckpt",
        device="cuda" if config["use_gpu"] else "cpu",
        cmd_copy_inputs=get_cmd_copy_inputs_synthseg,
    output:
        synthseg_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="synthseg",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids_log(
            "run_inference_synthseg",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    shadow:
        "minimal"
    threads: 8
    resources:
        gpus=1 if config["use_gpu"] else 0,
        mem_mb=16000,
        time=15 if config["use_gpu"] else 60,
    group:
        "subj"
    conda:
        "../envs/synthseg.yaml"
    shell:
        # Create temp model directory
        "mkdir -p {params.model_dir} && "
        # Extract model tar
        "tar -xf {input.model_tar} -C {params.model_dir} && "
        # Copy input
        "{params.cmd_copy_inputs} && "
        # Run SynthSeg inference
        "python {workflow.basedir}/scripts/seg_synthseg.py "
        "input.nii.gz "
        "{params.checkpoint_path} "
        "--output {output.synthseg_seg} "
        "--device {params.device} "
        "&> {log}"
