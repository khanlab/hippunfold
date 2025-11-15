"""nnUNet v1 inference rules"""
import re
from pathlib import Path


def get_nnunet_input(wildcards):
    T1w_nii = bids(
        root=root,
        datatype="anat",
        space="corobl",
        desc="preproc",
        hemi="{hemi}",
        suffix="T1w.nii.gz",
        **inputs.subj_wildcards,
    )
    T2w_nii = bids(
        root=root,
        datatype="anat",
        space="corobl",
        desc="preproc",
        hemi="{hemi}",
        suffix="T2w.nii.gz",
        **inputs.subj_wildcards,
    )
    if (config["modality"] == "T1w" or config["modality"] == "T2w") and config[
        "force_nnunet_model"
    ] == "T1T2w":
        return (T1w_nii, T2w_nii)

    elif config["modality"] == "T2w":
        return T2w_nii
    elif config["modality"] == "T1w":
        return T1w_nii
    elif config["modality"] == "hippb500":
        return bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            space="corobl",
            suffix="b500.nii.gz",
            **inputs.subj_wildcards,
        )
    else:
        raise ValueError("modality not supported for nnunet!")


def parse_task_from_tar(wildcards, input):
    match = re.search(r"Task[0-9]{3}_[\w]+", str(input.model_tar))
    if match:
        task = match.group(0)
    else:
        raise ValueError("cannot parse Task from model tar")
    return task


def parse_chkpnt_from_tar(wildcards, input):
    match = re.search(r"^.*\.(\w+)\.tar", str(input.model_tar))
    if match:
        chkpnt = match.group(1)
    else:
        raise ValueError("cannot parse chkpnt from model tar")
    return chkpnt


def parse_trainer_from_tar(wildcards, input):
    match = re.search(r"^.*\.(\w+)\..*.tar", str(input.model_tar))
    if match:
        trainer = match.group(1)
    else:
        raise ValueError("cannot parse trainer from model tar")
    return trainer


def get_cmd_copy_inputs(wildcards, input):
    in_img = input.in_img
    if isinstance(in_img, str):
        # we have one input image
        return f"cp {in_img} tempimg/temp_0000.nii.gz"
    else:
        cmd = []
        # we have multiple input images
        for i, img in enumerate(input.in_img):
            cmd.append(f"cp {img} tempimg/temp_{i:04d}.nii.gz")
        return " && ".join(cmd)


rule run_inference_v1:
    """nnUNet v1 inference with tar extraction in shadow directory"""
    input:
        in_img=get_nnunet_input,
        model_tar=get_model_tar(),
    params:
        cmd_copy_inputs=get_cmd_copy_inputs,
        temp_lbl="templbl/temp.nii.gz",
        model_dir="tempmodel",
        in_folder="tempimg",
        out_folder="templbl",
        task=parse_task_from_tar,
        chkpnt=parse_chkpnt_from_tar,
        trainer=parse_trainer_from_tar,
        tta="" if config["nnunet_enable_tta"] else "--disable_tta",
    output:
        nnunet_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="nnunetv1",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids_log(
            "run_inference_v1",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    shadow:
        "minimal"
    threads: 16
    resources:
        gpus=1 if config["use_gpu"] else 0,
        mem_mb=48000,
        time=30 if config["use_gpu"] else 120,
    group:
        "subj"
    conda:
        "../envs/nnunet.yaml"
    shell:
        "mkdir -p {params.model_dir} {params.in_folder} {params.out_folder} && "
        "{params.cmd_copy_inputs} && "
        "tar -xf {input.model_tar} -C {params.model_dir} && "
        "export RESULTS_FOLDER={params.model_dir} && "
        "export nnUNet_n_proc_DA={threads} && "
        "nnUNet_predict -i {params.in_folder} -o {params.out_folder} -t {params.task} -chk {params.chkpnt} -tr {params.trainer} {params.tta} &> {log} && "
        "cp {params.temp_lbl} {output.nnunet_seg}"
