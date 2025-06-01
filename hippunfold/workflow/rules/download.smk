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
    script:
        "../scripts/download.py"


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
        "../envs/c3d.yaml"
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
        "../envs/c3d.yaml"
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
        template_coords=temp(
            bids(
                root=root,
                datatype="coords",
                **inputs.subj_wildcards,
                dir="{dir}",
                label="{label}",
                suffix="coords.nii.gz",
                desc="init",
                space="template",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
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
        template_anat=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                hemi="{hemi}",
                suffix="{modality}.nii.gz".format(
                    modality=get_modality_suffix(config["modality"])
                ),
            ),
        ),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
    shell:
        "{params.copy_or_flip_cmd} {output.template_anat}"
