def get_model_tar():
    if config["force_nnunet_model"]:
        model_name = config["force_nnunet_model"]
    else:
        model_name = config["modality"]

    local_tar = config["nnunet_model"].get(model_name, None)
    if local_tar == None:
        print(f"ERROR: {model_name} does not exist in nnunet_model in the config file")

    return os.path.abspath(os.path.join(download_dir, local_tar.split("/")[-1]))


rule download_model:
    params:
        url=config["nnunet_model"][config["force_nnunet_model"]]
        if config["force_nnunet_model"]
        else config["nnunet_model"][config["modality"]],
    output:
        model_tar=get_model_tar(),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_tar}"


rule download_atlas:
    params:
        url=lambda wildcards: config["atlas_files_osf"][wildcards.atlas],
    output:
        model_zip=os.path.join(download_dir,"{atlas}"+'.zip')
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


def atlas_outs():
    outs = []
    for a in config["atlas"]:
        for fn in config["atlas_files"][a]:
            outs.append(os.path.join(download_dir,config["atlas_files"][a][fn]))
    return outs

rule unzip_atlas:
    input:
        model_zip=lambda wildcards: os.path.join(download_dir,"{wildcards.atlas}"+'.zip'),
        atlas=config["atlas"],
    params:
        dir=download_dir,
    output:
        atlas_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"


rule download_template:
    params:
        url=config["template_files_osf"][config["template"]],
    output:
        model_zip=os.path.join(download_dir,config["template"]+'.zip')
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


rule unzip_template:
    input:
        model_zip=os.path.join(download_dir,config["template"]+'.zip'),
    params:
        dir=download_dir,
    output:
        os.path.join(download_dir,config["template_files"][config["template"]]["T1w"])
    shell:
        "unzip {input.model_zip} -d {params.dir}"
