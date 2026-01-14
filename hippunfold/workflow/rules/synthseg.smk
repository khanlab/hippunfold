"""SynthSeg v1 inference rules"""
import re
from pathlib import Path


def get_synthseg_input(wildcards):
    """Get the appropriate input image for SynthSeg inference"""
    
    # SynthSeg typically works with T1w or T2w images
    if config["modality"] == "T2w":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="T2w.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "T1w":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="preproc", 
            hemi="{hemi}",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "hippb500":
        return bids(
            root=work,
            datatype="dwi",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="b500.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "dsegtissue":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="cropped",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    else:
        raise ValueError(f"Unsupported modality for SynthSeg: {config['modality']}")


def get_synthseg_input_for_flip(wildcards):
    """Get the appropriate input image for flipping before SynthSeg inference (left hemi only)"""
    
    if config["modality"] == "T2w":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="T2w.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "T1w":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="preproc", 
            hemi="{hemi}",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "hippb500":
        return bids(
            root=work,
            datatype="dwi",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
            suffix="b500.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    elif config["modality"] == "dsegtissue":
        return bids(
            root=work,
            datatype="anat",
            space="corobl",
            desc="cropped",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ).format(**wildcards)
    else:
        raise ValueError(f"Unsupported modality for SynthSeg: {config['modality']}")


# Rule to flip left hemisphere input before SynthSeg inference
rule flip_left_hemi_for_synthseg:
    """Flip left hemisphere image along x-axis before SynthSeg inference.
    
    SynthSeg models are typically trained on right hemisphere data,
    so we flip left hemisphere to match the expected orientation.
    """
    input:
        in_img=get_synthseg_input_for_flip,
    output:
        flipped_img=temp(
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="flipped.nii.gz",
                desc="synthseginput",
                space="corobl",
                hemi="{hemi,L}",
            )
        ),
    log:
        bids(
            root="logs",
            desc="flip_left_hemi_for_synthseg",
            **config["subj_wildcards"],
            hemi="{hemi}",
            suffix="log.txt"
        ),
    group:
        "subj"
#    conda:
#        "../envs/c3d.yaml"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.in_img} -flip x -o {output.flipped_img} &> {log}"


def get_cmd_copy_inputs_synthseg(wildcards, input):
    """Copy input images with SynthSeg naming convention"""
    in_img = input.in_img
    input_nii = "input.nii.gz"
    
    cmd = f"cp {in_img} {input_nii}"
    
    return cmd


def get_synthseg_inference_input(wildcards):
    """Get input for SynthSeg inference - flipped for left hemi, original for right hemi"""
    if wildcards.hemi == "L":
        # Left hemisphere: use flipped image
        return bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="flipped.nii.gz",
            desc="synthseginput",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    else:
        # Right hemisphere: use original image
        return get_synthseg_input(wildcards)


rule run_inference_synthseg:
    """SynthSeg inference with checkpoint extraction in shadow directory.
    
    For left hemisphere: input is pre-flipped, output will be unflipped in next rule.
    For right hemisphere: input and output are used as-is.
    """
    input:
        in_img=get_synthseg_inference_input,
        model_tar=get_model_tar(),
        script=str((Path(workflow.basedir) / "scripts" / "seg_synthseg.py").resolve())
    params:
        model_dir="tempmodel",
        checkpoint_path="tempmodel/synthseg/last.ckpt",
        cmd_copy_inputs=get_cmd_copy_inputs_synthseg,
        device="cuda" if config["use_gpu"] else "cpu",
    output:
        synthseg_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="dseg.nii.gz",
                desc="synthsegraw",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids(
            root="logs",
            desc="run_inference_synthseg",
            **config["subj_wildcards"],
            hemi="{hemi}",
            suffix="log.txt"
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
    shell:
        # Create temp model directory
        "mkdir -p {params.model_dir} && "
        # Extract model tar
        "tar -xf {input.model_tar} -C {params.model_dir} && "
        # Copy input
        "{params.cmd_copy_inputs} && "
        # Run SynthSeg inference
        "python3 {input.script} "
        "input.nii.gz "
        "{params.checkpoint_path} "
        "--output {output.synthseg_seg} "
        "--device {params.device} "
        "&> {log}"


# Rule to unflip left hemisphere prediction after SynthSeg inference
rule unflip_left_hemi_synthseg_output:
    """Unflip left hemisphere SynthSeg output along x-axis.
    
    After SynthSeg inference on flipped left hemisphere,
    we flip the prediction back to original orientation.
    """
    input:
        seg=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="synthsegraw",
            space="corobl",
            hemi="{hemi,L}",
        ),
    output:
        unflipped_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="dseg.nii.gz",
                desc="synthseg",
                space="corobl",
                hemi="{hemi,L}",
            )
        ),
    log:
        bids(
            root="logs",
            desc="unflip_left_hemi_synthseg_output",
            **config["subj_wildcards"],
            hemi="{hemi}",
            suffix="log.txt"
        ),
    group:
        "subj"
#    conda:
#        "../envs/c3d.yaml"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.seg} -flip x -o {output.unflipped_seg} &> {log}"


# Rule to pass through right hemisphere without flipping
rule passthrough_right_hemi_synthseg_output:
    """Pass through right hemisphere SynthSeg output without modification.
    
    Right hemisphere doesn't need flipping, so we just rename the file.
    """
    input:
        seg=bids(
            root=root,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="synthsegraw",
            space="corobl",
            hemi="{hemi,R}",
        ),
    output:
        seg=temp(
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="dseg.nii.gz",
                desc="synthseg",
                space="corobl",
                hemi="{hemi,R}",
            )
        ),
    group:
        "subj"
    shell:
        "cp {input.seg} {output.seg}"
