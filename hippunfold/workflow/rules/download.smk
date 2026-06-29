# populate the HIPPUNFOLD_CACHE_DIR folder as needed
from lib import utils as utils

download_dir = utils.get_download_dir()


localrules:
    import_template_dseg,
    import_template_dseg_dentate,
    import_template_coords,
    import_template_anat,
    import_template_anat_crop,
    cp_atlas_surf_gii,
    cp_atlas_metric_gii,
    download_extract_template,
    download_surf_template_atlas,


rule download_extract_template:
    """Note: OSF urls don't seem to be supported with snakemake storage plugin"""
    params:
        url=lambda wildcards: config["resource_urls"]["template"][wildcards.template],
    output:
        unzip_dir=directory(Path(download_dir) / "template" / "{template}"),
    group:
        "download"
    script:
        "../scripts/download.py"


rule download_surf_template_atlas:
    input:
        zip_file=lambda wildcards: storage(
            config["resource_urls"]["atlas"][wildcards.atlas]
        ),
    output:
        unzip_dir=temp(directory(Path(download_dir) / "atlases_dl" / "tpl-{atlas}")),
    wildcard_constraints:
        atlas="|".join(config["builtin_atlases"]),
    shell:
        "unzip {input} -d {output}"


rule cp_atlas_surf_gii:
    input:
        unzip_dir=Path(download_dir) / "atlases_dl" / "tpl-{atlas}",
    params:
        path=lambda wildcards, input: bids_atlas(
            root=Path(input.unzip_dir).parent,
            template=wildcards.atlas,
            hemi=wildcards.hemi,
            label=wildcards.label,
            den=wildcards.density,
            space=wildcards.space,
            suffix=f"{wildcards.surf_name}.surf.gii",
        ),
    output:
        atlas_file=bids_atlas(
            root=get_atlas_dir(),
            template="{atlas}",
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="{space}",
            suffix="{surf_name}.surf.gii",
        ),
    group:
        "atlas_and_qc"
    shell:
        "cp {params.path} {output}"


rule cp_atlas_metric_gii:
    input:
        unzip_dir=Path(download_dir) / "atlases_dl" / "tpl-{atlas}",
    params:
        path=lambda wildcards, input: bids_atlas(
            root=Path(input.unzip_dir).parent,
            template=wildcards.atlas,
            hemi=wildcards.hemi,
            label=wildcards.label,
            den=wildcards.density,
            suffix=f"{wildcards.metricname}.{wildcards.metrictype}.gii",
        ),
    output:
        atlas_file=bids_atlas(
            root=get_atlas_dir(),
            template="{atlas}",
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            suffix="{metricname}.{metrictype}.gii",
        ),
    group:
        "atlas_and_qc"
    shell:
        "cp {params.path} {output}"


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
    conda:
        "../envs/c3d.yaml"
    group:
        "download"
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
    conda:
        "../envs/c3d.yaml"
    group:
        "download"
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
    conda:
        "../envs/c3d.yaml"
    shell:
        "{params.copy_or_flip_cmd} {output.template_anat}"


rule import_template_anat_crop:  # used only in templateseg workflow
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_anat=lambda wildcards: Path(download_dir)
        / "template"
        / config["inject_template"]
        / config["template_files"][config["inject_template"]]["crop_ref"].format(
            **wildcards
        ),
        copy_or_flip_cmd=lambda wildcards: copy_or_flip(
            wildcards,
            Path(download_dir)
            / "template"
            / config["inject_template"]
            / config["template_files"][config["inject_template"]]["crop_ref"].format(
                **wildcards
            ),
        ),
    output:
        template_anat=temp(
            bids(
                root=root,
                datatype="anat",
                desc="template",
                space="corobl",
                **inputs.subj_wildcards,
                hemi="{hemi}",
                suffix="{modality}.nii.gz".format(
                    modality=get_modality_suffix(config["modality"])
                ),
            ),
        ),
    conda:
        "../envs/c3d.yaml"
    shell:
        "{params.copy_or_flip_cmd} {output.template_anat}"
