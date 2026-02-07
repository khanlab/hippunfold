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
    model_cfg = config["nnunet_models"].get(model_name, None)
    if model_cfg is None:
        raise KeyError(
            f"Model '{model_name}' does not exist in 'nnunet_models' in the config file"
        )

    url = model_cfg.get("url")
    if not url:
        raise ValueError(
            f"Model '{model_name}' in 'nnunet_models' is missing a valid 'url' entry"
        )
    tarfile = str((Path(download_dir) / "model" / Path(url).name).absolute())
    return tarfile


def get_model_dir():
    return get_model_tar().removesuffix(".gz").removesuffix(".tar")


rule download_nnunet_model:
    input:
        url=storage(model_dict["url"]),
    output:
        model_tar=get_model_tar(),
    container:
        config["singularity"]["autotop"]
    shell:
        "mkdir -p {params.model_dir} && wget https://{params.url} -O {output}"


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


rule unpack_nnunet_model:
    """ Unpack nnunet model tar to temp folder to check contents"""
    input:
        tar=get_model_tar(),
    params:
        tar_opts=lambda wildcards, input: "-xzf" if input.tar[-2:] == "gz" else "-xf",
    output:
        directory(get_model_dir()),
    shell:
        "mkdir -p {output} && tar {params.tar_opts} {input} -C {output}"


if model_dict["arch_version"] == "nnunet_v1":

    rule run_inference_nnunet_v1:
        """This rule uses either GPU or CPU .
        It also runs in an isolated folder (shadow), with symlinks to inputs in that folder, copying over outputs once complete, so temp files are not retained
        """
        input:
            in_img=get_nnunet_input,
            model_dir=get_model_dir(),
        params:
            cmd_copy_inputs=get_cmd_copy_inputs,
            temp_lbl="templbl/temp.nii.gz",
            in_folder="tempimg",
            out_folder="templbl",
            task=model_dict["task"],
            chkpnt=model_dict["checkpoint"],
            trainer=model_dict["trainer"],
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
        conda:
            "../envs/nnunet.yaml"
        shell:
            "mkdir -p {params.in_folder} {params.out_folder} && "
            "{params.cmd_copy_inputs} && "
            "export RESULTS_FOLDER={input.model_dir} && "
            "export nnUNet_n_proc_DA={threads} && "
            "nnUNet_predict -i {params.in_folder} -o {params.out_folder} -t {params.task} -chk {params.chkpnt} -tr {params.trainer} {params.tta} &> {log} && "
            "cp {params.temp_lbl} {output.nnunet_seg}"

elif model_dict["arch_version"] == "nnunet_v2":

    rule run_inference_nnunet_v2:
        """nnUNet v2 inference with tar extraction in shadow directory"""
        input:
            in_img=get_nnunet_input,
            model_tar=get_model_tar(),
        params:
            cmd_copy_inputs=get_cmd_copy_inputs,
            tar_opts=lambda wildcards, input: (
                "-xzf" if input.model_tar[-2:] == "gz" else "-xf"
            ),
            temp_lbl="templbl/temp.nii.gz",
            model_dir="tempmodel",
            in_folder="tempimg",
            out_folder="templbl",
            dataset_id=model_dict["dataset_id"],
            configuration=model_dict["configuration"],
            trainer=model_dict["trainer"],
            plans=model_dict["plans"],
            chkpnt=model_dict["checkpoint"],
            tta="" if config["nnunet_enable_tta"] else "--disable_tta",
            device="cuda" if config["use_gpu"] else "cpu",
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
                "run_inference_nnunet_v2",
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
            "../envs/nnunetv2.yaml"
        shell:
            "mkdir -p {params.model_dir} {params.in_folder} {params.out_folder} && "

            "{params.cmd_copy_inputs} && "

            "tar {params.tar_opts} {input.model_tar} -C {params.model_dir} && "

            "export nnUNet_results={params.model_dir}/nnunet_v2 && "
            "export nnUNet_raw={params.model_dir}/nnunet_v2/nnUNet_raw && "
            "export nnUNet_preprocessed={params.model_dir}/nnunet_v2/nnUNet_preprocessed && "
            "export nnUNet_n_proc_DA={threads} && "

            "FOLDS=$(find {params.model_dir}/nnunet_v2/Dataset{params.dataset_id}_*/nnUNetTrainer* -maxdepth 1 -type d -name 'fold_*' | sed 's/.*fold_//' | sort -n | tr '\\n' ' ' | sed 's/ $//' ) && "

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
            "cp {params.temp_lbl} {output.nnunet_seg}"

elif model_dict["arch_version"] == "synthseg_v2":

    rule flip_synthseg_input:
        input:
            nii=bids(
                root=root,
                datatype="{datatype}",
                suffix="{suffix}.nii.gz",
                desc="preproc",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            ),
        output:
            nii=temp(
                bids(
                    root=root,
                    datatype="{datatype}",
                    suffix="{suffix}.nii.gz",
                    desc="preproc",
                    space="corobl",
                    hemi="{hemi}flip",
                    **inputs.subj_wildcards,
                )
            ),
        conda:
            "../envs/c3d.yaml"
        group:
            "subj"
        shell:
            "c3d {input} -flip x {output}"

    rule run_inference_synthseg_v2:
        """SynthSeg inference with checkpoint extraction in shadow directory.

        For left hemisphere: input is pre-flipped, output will be unflipped in next rule.
        For right hemisphere: input and output are used as-is.
        """
        input:
            in_img=get_nnunet_input,
            model_tar=get_model_tar(),
        params:
            model_dir="tempmodel",
            checkpoint_path="tempmodel/synthseg/{chkpt}".format(
                chkpt=model_dict["checkpoint"]
            ),
            device="cuda" if config["use_gpu"] else "cpu",
        output:
            synthseg_seg=temp(
                bids(
                    root=root,
                    datatype="anat",
                    suffix="dseg.nii.gz",
                    desc="nnunet",
                    space="corobl",
                    hemi="{hemi,Lflip|R}",
                    **inputs.subj_wildcards,
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

            "tar -xf {input.model_tar} -C {params.model_dir} && "

            "python {workflow.basedir}/scripts/seg_synthseg.py "
            "{input.in_img} "
            "{params.checkpoint_path} "
            "--output {output.synthseg_seg} "
            "--device {params.device} "
            "&> {log}"
            # Extract model tar
            # Run SynthSeg inference

    rule unflip_synthseg_output:
        input:
            nii=bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="nnunet",
                space="corobl",
                hemi="{hemi}flip",
                **inputs.subj_wildcards,
            ),
        output:
            nii=temp(
                bids(
                    root=root,
                    datatype="anat",
                    suffix="dseg.nii.gz",
                    desc="nnunet",
                    space="corobl",
                    hemi="{hemi,L}",
                    **inputs.subj_wildcards,
                )
            ),
        conda:
            "../envs/c3d.yaml"
        group:
            "subj"
        shell:
            "c3d {input} -flip x {output}"


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
