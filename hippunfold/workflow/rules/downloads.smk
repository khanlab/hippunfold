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

rule download_template:
    params:
        url=lambda wildcards: config["template_files_osf"][wildcards.atlas],
    output:
        model_zip=os.path.join(download_dir,"{template}"+'.zip')
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

def template_outs():
    outs = []
    for fn in config["template_files"][config["template"]]:
        outs.append(os.path.join(download_dir,config["template_files"][config["template"]][fn]))
    return outs


rule unzip_template:
    input:
        model_zip=os.path.join(download_dir,"{template}"+'.zip'),
    params:
        dir=download_dir,
    output:
        template_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"

rule unzip_atlas:
    input:
        model_zip=os.path.join(download_dir,"{atlas}"+'.zip'),
    params:
        dir=download_dir,
    output:
        atlas_outs()
    shell:
        "unzip {input.model_zip} -d {params.dir}"


rule unzip_template_shape:
    input:
        model_zip=os.path.join(download_dir,config["inject_template"]+'.zip'),
    params:
        dir=download_dir,
    output:
        template_seg=os.path.join(
            download_dir,
            config["template_files"][config["inject_template"]]["dseg"],
        )
    shell:
        "unzip {input.model_zip} -d {params.dir}"
    
