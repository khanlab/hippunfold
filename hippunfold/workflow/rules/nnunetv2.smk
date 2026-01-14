"""nnUNet v2 inference rules"""
import re
from pathlib import Path


def parse_dataset_from_tar_v2(wildcards, input):
    """Parse Dataset### from nnUNet v2 model tar filename"""
    match = re.search(r"Dataset([0-9]{3})_[\w]+", str(input.model_tar))
    if match:
        dataset_id = match.group(1)
        return dataset_id
    else:
        raise ValueError("Cannot parse Dataset from model tar for nnUNet v2")


def parse_configuration_from_tar_v2(wildcards, input):
    """Parse configuration (e.g., 3d_fullres) from nnUNet v2 model tar filename"""
    # v2 structure: Dataset###_name/trainer__plans__config/
    # Look for patterns like __3d_fullres, __2d, __3d_lowres in tar name
    match = re.search(r"__([23]d(?:_\w+)?)", str(input.model_tar))
    if match:
        return match.group(1)
    # Default to 3d_fullres if not found
    return "3d_fullres"


def parse_trainer_from_tar_v2(wildcards, input):
    """Parse trainer from nnUNet v2 model tar filename"""
    # v2 tar name pattern: Dataset###_name.trainer__plans__config.tar.gz
    match = re.search(r"Dataset[0-9]{3}_[\w]+\.([\w]+)__", str(input.model_tar))
    if match:
        return match.group(1)
    return "nnUNetTrainer"


def parse_plans_from_tar_v2(wildcards, input):
    """Parse plans identifier from nnUNet v2 model tar filename"""
    # v2 tar name pattern: Dataset###_name.trainer__plans__config.tar.gz
    match = re.search(r"__([\w]+)__[23]d", str(input.model_tar))
    if match:
        return match.group(1)
    return "nnUNetPlans"


def get_checkpoint_name_v2(wildcards, input):
    """
    Parse checkpoint name from v2 model tar filename.
    nnUNet v2 will automatically use whatever checkpoint is available if not found.
    """
    # Try to parse checkpoint type from tar filename
    # Pattern: Dataset###_name.trainer__plans__config.checkpoint_XXX.tar.gz
    match = re.search(r"\.(checkpoint_\w+)\.tar", str(input.model_tar))
    if match:
        return match.group(1) + ".pth"
    
    # If no checkpoint specified in filename, return "checkpoint_best.pth" as default
    # nnUNet v2 will auto-detect available checkpoints if this doesn't exist
    return "checkpoint_best.pth"


def get_cmd_copy_inputs_v2(wildcards, input):
    """Copy input images with nnUNet v2 naming convention"""
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


rule run_inference_v2:
    """nnUNet v2 inference with tar extraction in shadow directory"""
    input:
        in_img=get_nnunet_input,  # Use same input function from nnunet.smk
        model_tar=get_model_tar(),
    params:
        cmd_copy_inputs=get_cmd_copy_inputs_v2,
        temp_lbl="templbl/temp.nii.gz",
        model_dir="tempmodel",
        in_folder="tempimg",
        out_folder="templbl",
        dataset_id=parse_dataset_from_tar_v2,
        configuration=parse_configuration_from_tar_v2,
        trainer=parse_trainer_from_tar_v2,
        plans=parse_plans_from_tar_v2,
        chkpnt=get_checkpoint_name_v2,
        tta="" if config["nnunet_enable_tta"] else "--disable_tta",
        device="cuda" if config["use_gpu"] else "cpu",
    output:
        nnunet_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="dseg.nii.gz",
                desc="nnunetv2",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids(
            root="logs",
            desc="run_inference_v2",
            **config["subj_wildcards"],
            hemi="{hemi}",
            suffix="log.txt"
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
    container:
        config["singularity"]["nnunetv2"]
    shell:
        # Create temp folders
        "mkdir -p {params.model_dir} {params.in_folder} {params.out_folder} && "
        # Copy input images
        "{params.cmd_copy_inputs} && "
        # Extract model tar
        "tar -xf {input.model_tar} -C {params.model_dir} && "
        # Set nnUNet v2 environment variables - handle nested nnunet directory structure
        "export nnUNet_results={params.model_dir}/nnunet_v2 && "
        "export nnUNet_raw={params.model_dir}/nnunet_v2/nnUNet_raw && "
        "export nnUNet_preprocessed={params.model_dir}/nnunet_v2/nnUNet_preprocessed && "
        "export nnUNet_n_proc_DA={threads} && "
        # Detect available folds
        "FOLDS=$(find {params.model_dir}/nnunet_v2/Dataset{params.dataset_id}_*/nnUNetTrainer* -maxdepth 1 -type d -name 'fold_*' | sed 's/.*fold_//' | sort -n | tr '\\n' ' ' | sed 's/ $//' ) && "
        # Run nnUNet v2 prediction with explicit device specification
        "nnUNetv2_predict "
        "-i {params.in_folder} "
        "-o {params.out_folder} "
        "-d {params.dataset_id} "
        "-c {params.configuration} "
        "-tr {params.trainer} "
        "-p {params.plans} "
        "-f $FOLDS "
        "-chk {params.chkpnt} "
        "-device {params.device} "
        "{params.tta} "
        "&> {log} && "
        # Copy output
        "cp {params.temp_lbl} {output.nnunet_seg}"
