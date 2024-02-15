# populate the HIPPUNFOLD_CACHE_DIR folder as needed

download_dir = get_download_dir()


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
