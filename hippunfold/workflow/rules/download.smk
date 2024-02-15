# populate the HIPPUNFOLD_CACHE_DIR folder as needed

download_dir = get_download_dir()


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


rule download_extract_atlas_or_template:
    params:
        url=lambda wildcards: config["resource_urls"][wildcards.resource_type][
            wildcards.atlas
        ],
    output:
        unzip_dir=directory(
            Path(download_dir) / "{resource_type,atlas|template}" / "{atlas}"
        ),
    container:
        config["singularity"]["autotop"]
    shadow:
        "minimal"
    shell:
        "wget https://{params.url} -O temp.zip && "
        " unzip -d {output.unzip_dir} temp.zip"
