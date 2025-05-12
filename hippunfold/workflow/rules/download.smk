# populate the HIPPUNFOLD_CACHE_DIR folder as needed
from lib import utils as utils

download_dir = utils.get_download_dir()


rule download_extract_template:
    params:
        url=lambda wildcards: config["resource_urls"]["template"][wildcards.template],
    output:
        unzip_dir=directory(Path(download_dir) / "template" / "{template}"),
    shadow:
        "minimal"
    conda:
        conda_env("curl")
    shell:
        "curl -L 'https://{params.url}' -o temp.zip && "
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
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="hipptissue",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "{params.copy_or_flip_cmd} {output.template_seg}"


rule import_template_layers:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_seg=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]]["layers"].format(
            **wildcards
        ),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]]["layers"].format(
                **wildcards
            ),
        ),
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="hipplayers",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "{params.copy_or_flip_cmd} {output.template_seg}"


rule import_template_dseg_dentate:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_seg=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]]["dseg_dentate"].format(
            **wildcards
        ),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]][
                "dseg_dentate"
            ].format(**wildcards),
        ),
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="dentatetissue",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "{params.copy_or_flip_cmd} {output.template_seg}"

