# populate the HIPPUNFOLD_CACHE_DIR folder as needed

download_dir = get_download_dir()


rule download_extract_atlas:
    params:
        url=lambda wildcards: config["resource_urls"]["atlas"][wildcards.atlas],
    output:
        unzip_dir=directory(Path(download_dir) / "atlas" / "{atlas}"),
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"
    shell:
        "wget 'https://{params.url}' -O temp.zip && "
        " unzip -d {output.unzip_dir} temp.zip"


rule download_extract_template:
    params:
        url=lambda wildcards: config["resource_urls"]["template"][wildcards.template],
    output:
        unzip_dir=directory(Path(download_dir) / "template" / "{template}"),
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"
    shell:
        "wget 'https://{params.url}' -O temp.zip && "
        " unzip -d {output.unzip_dir} temp.zip"


## unpack template
# this is used in both shape_inject.smk and templateseg.smk, so instead of having redundant rules I will simply apply them hemisphere
# Template-based segmentation supports templates that have only a single hemisphere
# by flipping it


def copy_or_flip(wildcards, file_to_process):
    if (
        wildcards.hemi
        in config["template_files"][config["inject_template"]]["hemi_wildcards"]
    ):
        cmd = f"cp {file_to_process}"
    else:
        cmd = f"c3d {file_to_process} -flip x -o"
    return cmd


rule import_template_dseg:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_seg=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]]["dseg"].format(
            **wildcards
        ),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]]["dseg"].format(
                **wildcards
            ),
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "{params.copy_or_flip_cmd} {output.template_seg}"


rule import_template_coords:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_coords=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]]["coords"].format(
            **wildcards
        ),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]]["coords"].format(
                **wildcards
            ),
        ),
    output:
        template_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi="{hemi}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "{params.copy_or_flip_cmd} {output.template_coords}"


rule import_template_anat:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_anat=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]][
            get_modality_suffix(config["modality"])
        ].format(**wildcards),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]][
                get_modality_suffix(config["modality"])
            ].format(**wildcards),
        ),
    output:
        template_anat=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "{params.copy_or_flip_cmd} {output.template_anat}"


def get_model_tar():
    if config["force_nnunet_model"]:
        model_name = config["force_nnunet_model"]
    else:
        model_name = config["modality"]

    local_tar = config["resource_urls"]["nnunet_model"].get(model_name, None)
    if local_tar == None:
        print(f"ERROR: {model_name} does not exist in nnunet_model in the config file")

    return (Path(download_dir) / "model" / Path(local_tar).name).absolute()


rule download_extract_nnunet_model:
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
        "mkdir -p {params.model_dir} && wget https://{params.url} -O {output} && "
        "tar -xf {output} -C {params.model_dir}"
