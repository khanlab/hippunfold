import re


def get_model_tar():
    model_name = config["nnunet_model"]

    local_tar = config["resource_urls"]["nnunet_model"].get(model_name, None)
    if local_tar == None:
        print(f"ERROR: {model_name} does not exist in nnunet_model in the config file")

    return (Path(download_dir) / "model" / Path(local_tar).name).absolute()


rule download_nnunet_model:
    params:
        url=(
            config["resource_urls"]["nnunet_model"][config["nnunet_model"]]
        ),
        model_dir=Path(download_dir) / "model",
    output:
        model_tar=get_model_tar(),
    conda:
        conda_env("curl")
    shell:
        "mkdir -p {params.model_dir} && curl -L https://{params.url} -o {output}"


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


rule run_inference:
    """ This rule uses either GPU or CPU .
    It also runs in an isolated folder (shadow), with symlinks to inputs in that folder, copying over outputs once complete, so temp files are not retained"""
    input:
        in_img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix="preproc.nii.gz",
            **inputs.subj_wildcards,
        ),
        model_tar=get_model_tar(),
    params:
        cmd_copy_inputs=get_cmd_copy_inputs,
        temp_lbl="templbl/temp.nii.gz",
        model_dir="tempmodel",
        in_folder="tempimg",
        out_folder="templbl",
        tta="" if config["nnunet_enable_tta"] else "--disable_tta",
    output:
        nnunet_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="nnunet",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids_log(
            "run_inference",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    shadow:
        "minimal"
    threads: 16
    resources:
        gpus=1 if config["use_gpu"] else 0,
        mem_mb=16000,
        time=30 if config["use_gpu"] else 60,
    group:
        "subj"
    # conda:
    #     conda_env("nnunet")
    shell:
        #create temp folders
        #cp input image to temp folder
        #extract model
        #set nnunet env var to point to model
        #set threads
        # run inference
        #copy from temp output folder to final output
        "mkdir -p {params.model_dir} {params.in_folder} {params.out_folder} && "
        "{params.cmd_copy_inputs} && "
        "tar -xf {input.model_tar} -C {params.model_dir} && "
        "export nnUNet_results={params.model_dir} && "
        "export nnUNet_n_proc_DA={threads} && "
        "nnUNetv2_predict -i {params.in_folder} -o {params.out_folder} -d 001 -c 3d_fullres {params.tta} &> {log} && "
        "cp {params.temp_lbl} {output.nnunet_seg}"


rule qc_nnunet_dice:
    input:
        res_mask=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="nnunet",
                space="corobl",
                hemi="{hemi}",
            )
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
        conda_env("pyunfold")
    script:
        "../scripts/dice.py"
