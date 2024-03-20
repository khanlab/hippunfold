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
        "wget https://{params.url} -O temp.zip && "
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
        "wget https://{params.url} -O temp.zip && "
        " unzip -d {output.unzip_dir} temp.zip"
