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
        url=config["atlas_links_osf"][config["atlas"]],
    output:
        model_zip=os.path.join(download_dir, config["atlas"] + ".zip"),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


rule download_template:
    params:
        url=config["template_links_osf"][config["template"]],
    output:
        model_zip=os.path.join(download_dir, config["template"] + ".zip"),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


rule download_template_shape:
    params:
        url=config["template_links_osf"][config["inject_template"]],
    output:
        model_zip=os.path.join(download_dir, config["inject_template"] + ".zip"),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


def atlas_outs():
    outs = []
    for fn in config["atlas_files"][config["atlas"]]:
        for hemi in ["L", "R"]:
            outs.append(
                download_dir
                + "/"
                + expand(config["atlas_files"][config["atlas"]][fn], hemi=hemi)[0]
            )
    return list(dict.fromkeys(outs))


def template_outs():
    outs = []
    for fn in config["template_files"][config["template"]]:
        for hemi in ["L", "R"]:
            for dir in ["AP", "PD", "IO"]:
                outs.append(
                    download_dir
                    + "/"
                    + expand(
                        config["template_files"][config["template"]][fn],
                        hemi=hemi,
                        dir=dir,
                    )[0]
                )
    return list(dict.fromkeys(outs))


def template_shape_outs():
    outs = []
    for fn in config["template_files"][config["inject_template"]]:
        for hemi in ["L", "R"]:
            for dir in ["AP", "PD", "IO"]:
                for autotop in ["hipp", "dentate"]:
                    outs.append(
                        download_dir
                        + "/"
                        + expand(
                            config["template_files"][config["inject_template"]][fn],
                            hemi=hemi,
                            dir=dir,
                            autotop=autotop,
                        )[0]
                    )
    return list(dict.fromkeys(outs))


rule unzip_template:
    input:
        model_zip=os.path.join(download_dir, config["template"] + ".zip"),
    params:
        dir=os.path.join(download_dir, config["template"]),
    output:
        template_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"


rule unzip_atlas:
    input:
        model_zip=os.path.join(download_dir, config["atlas"] + ".zip"),
    params:
        dir=os.path.join(download_dir, config["atlas"]),
    output:
        atlas_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"


rule unzip_template_shape:
    input:
        model_zip=os.path.join(download_dir, config["inject_template"] + ".zip"),
    params:
        dir=os.path.join(download_dir, config["inject_template"]),
    output:
        template_shape_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"
