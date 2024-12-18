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


## unpack template
# this is used in both shape_inject.smk and templateseg.smk, so instead of having redundant rules I will simply apply them hemisphere
# Template-based segmentation supports templates that have only a single hemisphere
# by flipping it

hemi_constraints = []
if config["template"] in config["template_available_hemis"]:
    for hemi in config["hemi"]:
        if hemi in config["template_available_hemis"][config["template"]]:
            hemi_constraints.append(hemi)

hemi_constrained_wildcard = "{{hemi,{constraints}}}".format(
    constraints="|".join(hemi_constraints)
)


rule import_template_dseg:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_seg=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["inject_template"]]["dseg"].format(
            **wildcards
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi=hemi_constrained_wildcard,
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    shell:
        "cp {params.template_seg} {output.template_seg}"


rule flip_template_dseg:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi=hemi_constrained_wildcard,
            suffix="dseg.nii.gz"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="hipptissue",
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.template_seg} -flip x -o {output.nii} "


rule import_template_coords:
    input:
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        template_coords=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["coords"].format(**wildcards),
    output:
        template_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi=hemi_constrained_wildcard,
        ),
    group:
        "subj"
    shell:
        "cp {params.template_coords} {output.template_coords}"


rule flip_template_coords:
    input:
        template_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi=hemi_constrained_wildcard,
        ),
    output:
        template_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi="{hemi}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.template_coords} -flip x -o {output.template_coords} "


rule import_template_anat:
    input:
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        template_anat=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]][
            get_modality_suffix(config["modality"])
        ].format(**wildcards),
    output:
        template_anat=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            hemi=hemi_constrained_wildcard,
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
        ),
    group:
        "subj"
    shell:
        "cp {params.template_anat} {output.template_anat}"


rule flip_template_anat:
    input:
        template_anat=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            hemi=hemi_constrained_wildcard,
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
        ),
    output:
        template_anat=bids(
            root=work,
            datatype="anat",
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.template_anat} -flip x -o {output.template_anat} "
