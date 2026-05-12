import re


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
        model_tar=temp(get_model_tar()),
    shell:
        "cp {input} {output}"


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
                    desc="pred",
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
                    desc="pred",
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

elif model_dict["arch_version"] == "synthtopo":

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

    rule run_inference_synthtopo:
        """SynthTopo inference with checkpoint extraction in shadow directory.

        For left hemisphere: input is pre-flipped, output will be unflipped in next rule.
        For right hemisphere: input and output are used as-is.
        """
        input:
            in_img=get_nnunet_input,
            model_tar=get_model_tar(),
        params:
            model_dir="tempmodel",
            checkpoint_name=model_dict["checkpoint"],
            nb_classes=model_dict.get("nb_classes", 9),
            seg_nb_levels=model_dict.get("seg_nb_levels", 6),
            seg_features=" ".join(
                str(x) for x in model_dict.get("seg_features", [16, 24, 32, 48, 64, 96])
            ),
            device="cuda" if config["use_gpu"] else "cpu",
        output:
            synthseg_seg=temp(
                bids(
                    root=root,
                    datatype="anat",
                    suffix="dseg.nii.gz",
                    desc="pred",
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
            "../envs/synthtopo.yaml"
        shell:
            # Create temp model directory
            "mkdir -p {params.model_dir} && "

            # Extract model archive (tar.* or zip)
            "case \"{input.model_tar}\" in "
            "  *.zip) unzip -q \"{input.model_tar}\" -d \"{params.model_dir}\" ;; "
            "  *) tar -xf \"{input.model_tar}\" -C \"{params.model_dir}\" ;; "
            "esac && "

            # Locate checkpoint inside extracted tree
            "CHKPT=\"$(find \"{params.model_dir}\" -type f -name \"{params.checkpoint_name}\" -print -quit)\" && "
            "if [ -z \"$CHKPT\" ]; then echo \"ERROR: checkpoint not found: {params.checkpoint_name}\" 1>&2; exit 1; fi && "

            "python {workflow.basedir}/scripts/synthtopo.py "
            "{input.in_img} "
            "\"$CHKPT\" "
            "--output {output.synthseg_seg} "
            "--device {params.device} "
            "--nb-classes {params.nb_classes} "
            "--seg-nb-levels {params.seg_nb_levels} "
            "--seg-features {params.seg_features} "
            "&> {log}"
            # Extract model tar
            # Run SynthSeg inference

    rule unflip_synthseg_output:
        input:
            nii=bids(
                root=root,
                datatype="anat",
                suffix="dseg.nii.gz",
                desc="pred",
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
                    desc="pred",
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


def get_mask_crop_ref(wildcards, input):
    template_name = config.get("template")
    template_cfg = config.get("template_files", {}).get(template_name, {})
    mask_tpl = template_cfg.get("Mask_crop")
    if not mask_tpl:
        raise ValueError(
            f"Dice QC requires template_files.{template_name}.Mask_crop to be set; "
            "either add it to the config or omit the unetf3d dice QC output."
        )
    return Path(input.template_dir) / mask_tpl.format(**wildcards)


rule qc_nnunet_f3d:
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
            desc="pred",
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
        ref=get_mask_crop_ref,
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
