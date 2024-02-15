# create rules for downloading each atlas, since output files need to be specified
# without using wildcards
for atlas in config["atlas"]:

    rule:
        name:
            f"download_atlas_{atlas}"
        params:
            url=config["atlas_links_url"][atlas],
        output:
            model_zip=os.path.join(download_dir, atlas + ".zip"),
        container:
            config["singularity"]["autotop"]
        shell:
            "wget https://{params.url} -O {output.model_zip}"

    rule:
        name:
            f"unzip_download_atlas_{atlas}"
        input:
            model_zip=os.path.join(download_dir, atlas + ".zip"),
        params:
            dir=os.path.join(download_dir, atlas),
        output:
            [
                expand(Path(download_dir) / path, hemi=config["hemi"])
                for key, path in config["atlas_files"][atlas].items()
            ],
        shell:
            "unzip {input.model_zip} -d {params.dir}"


rule download_template:
    params:
        url=config["template_links_url"][config["template"]],
    output:
        model_zip=os.path.join(download_dir, config["template"] + ".zip"),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


rule download_template_shape:
    params:
        url=config["template_links_url"][config["inject_template"]],
    output:
        model_zip=os.path.join(download_dir, config["inject_template"] + ".zip"),
    container:
        config["singularity"]["autotop"]
    shell:
        "wget https://{params.url} -O {output.model_zip}"


def atlas_outs():
    outs = []
    for fn in config["atlas_files"][wildcards.atlas]:
        for hemi in ["L", "R"]:
            outs.append(
                download_dir
                + "/"
                + expand(config["atlas_files"][wildcards.atlas][fn], hemi=hemi)[0]
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


rule unzip_template_shape:
    input:
        model_zip=os.path.join(download_dir, config["inject_template"] + ".zip"),
    params:
        dir=os.path.join(download_dir, config["inject_template"]),
    output:
        template_shape_outs(),
    shell:
        "unzip {input.model_zip} -d {params.dir}"
