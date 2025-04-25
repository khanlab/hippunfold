import re


def get_nnunet_input(wildcards):
    T1w_nii = bids(
        root=work,
        datatype="anat",
        **config["subj_wildcards"],
        suffix="T1w.nii.gz",
        space="corobl",
        desc="preproc",
        hemi="{hemi}",
    )
    T2w_nii = bids(
        root=work,
        datatype="anat",
        **config["subj_wildcards"],
        suffix="T2w.nii.gz",
        space="corobl",
        desc="preproc",
        hemi="{hemi}",
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
            root=work,
            datatype="dwi",
            hemi="{hemi}",
            space="corobl",
            suffix="b500.nii.gz",
            **config["subj_wildcards"],
        )
    else:
        raise ValueError("modality not supported for nnunet!")


def get_model_tar():
    if config["force_nnunet_model"]:
        model_name = config["force_nnunet_model"]
    else:
        model_name = config["modality"]

    local_tar = config["resource_urls"]["nnunet_model"].get(model_name, None)
    if local_tar == None:
        print(f"ERROR: {model_name} does not exist in nnunet_model in the config file")

    return (Path(download_dir) / "model" / Path(local_tar).name).absolute()


rule download_nnunet_model:
    params:
        url=config["resource_urls"]["nnunet_model"][config["force_nnunet_model"]]
        if config["force_nnunet_model"]
        else config["resource_urls"]["nnunet_model"][config["modality"]],
        model_dir=Path(download_dir) / "model",
    output:
        model_tar=get_model_tar(),
    container:
        config["singularity"]["autotop"]
    shell:
        "mkdir -p {params.model_dir} && wget https://{params.url} -O {output}"


def parse_task_from_tar(wildcards, input):
    match = re.search("Task[0-9]{3}_[\w]+", input.model_tar)
    if match:
        task = match.group(0)
    else:
        raise ValueError("cannot parse Task from model tar")
    return task


def parse_chkpnt_from_tar(wildcards, input):
    match = re.search("^.*\.(\w+)\.tar", input.model_tar)
    if match:
        chkpnt = match.group(1)
    else:
        raise ValueError("cannot parse chkpnt from model tar")
    return chkpnt


def parse_trainer_from_tar(wildcards, input):
    match = re.search("^.*\.(\w+)\..*.tar", input.model_tar)
    if match:
        trainer = match.group(1)
    else:
        raise ValueError("cannot parse chkpnt from model tar")
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


rule run_inference:
    """ This rule uses either GPU or CPU .
    It also runs in an isolated folder (shadow), with symlinks to inputs in that folder, copying over outputs once complete, so temp files are not retained"""
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
        nnunet_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi,Lflip|R}"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="nnunet.txt",
            space="corobl",
            hemi="{hemi,Lflip|R}"
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
    container:
        config["singularity"]["autotop"]
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
        "export RESULTS_FOLDER={params.model_dir} && "
        "export nnUNet_n_proc_DA={threads} && "
        "nnUNet_predict -i {params.in_folder} -o {params.out_folder} -t {params.task} -chk {params.chkpnt} -tr {params.trainer} {params.tta} &> {log} && "
        "cp {params.temp_lbl} {output.nnunet_seg}"


rule unflip_nnunet_nii:
    """Unflip the Lflip nnunet seg"""
    input:
        nnunet_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}flip"
        ),
        unflip_ref=(
            bids(
                root=work,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="{modality}.nii.gz".format(modality=config["modality"]),
                space="corobl",
                desc="preproc",
                hemi="{hemi}",
            ),
        ),
    output:
        nnunet_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi,L}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nnunet_seg} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nnunet_seg} "


def get_f3d_ref(wildcards, input):
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
    input:
        img=(
            bids(
                root=work,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="{modality}.nii.gz".format(modality=config["modality"]),
                space="corobl",
                desc="preproc",
                hemi="{hemi}",
            ),
        ),
        seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=get_f3d_ref,
    output:
        cpp=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="cpp.nii.gz",
            desc="f3d",
            space="corobl",
            hemi="{hemi}"
        ),
        res=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{modality}.nii.gz".format(modality=config["modality"]),
            desc="f3d",
            space="template",
            hemi="{hemi}"
        ),
        res_mask=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="mask.nii.gz",
            desc="f3d",
            space="template",
            hemi="{hemi}"
        ),
    container:
        config["singularity"]["autotop"]
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="qcreg.txt",
            desc="f3d",
            space="corobl",
            hemi="{hemi}"
        ),
    group:
        "subj"
    shell:
        "reg_f3d -flo {input.img} -ref {params.ref} -res {output.res} -cpp {output.cpp} &> {log} && "
        "reg_resample -flo {input.seg} -cpp {output.cpp} -ref {params.ref} -res {output.res_mask} -inter 0 &> {log}"


rule qc_nnunet_dice:
    input:
        res_mask=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="mask.nii.gz",
            desc="f3d",
            space="template",
            hemi="{hemi}"
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
                **config["subj_wildcards"]
            ),
            caption="../report/nnunet_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/dice.py"
