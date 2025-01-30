def get_smoothing_opt(wildcards):
    """sets the smoothness of the greedy template shape injection deformation"""

    gradient_sigma = 1.732 * float(config["template_seg_smoothing_factor"])
    warp_sigma = 0.7071 * float(config["template_seg_smoothing_factor"])

    return f"-s {gradient_sigma}vox {warp_sigma}vox"


rule template_reg:
    input:
        fixed_img=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            space="corobl",
            desc="preproc",
            hemi="{hemi}"
        ),
        moving_img=bids(
            root=work,
            datatype="anat",
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        xfm_corobl=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
        general_opts="-d 3 -m NCC 2x2x2",
        smoothing_opts=get_smoothing_opt,
        iteration_opts="-n 100x50x10",  #default -n 100x100
    output:
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} "
        " {params.smoothing_opts} {params.iteration_opts} "
        " -i {input.fixed_img} {input.moving_img} -it {params.xfm_corobl} -o {output.warp}"


rule warp_template_dseg:
    input:
        ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
        template_dseg=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="hipptissue",
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    params:
        interp_opt="-ri LABEL 0.2vox",
    output:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {input.template_dseg} {output.inject_seg}  -r {input.warp}"


rule warp_template_coords:
    input:
        template_dir=Path(download_dir) / "template" / config["template"],
        ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
        template_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi="{hemi}"
        ),
    params:
        interp_opt="-ri NN",
    output:
        init_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {input.template_coords} {output.init_coords}  -r {input.warp}"


rule warp_template_anat:
    input:
        ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
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
    params:
        xfm_corobl=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
    output:
        warped=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            desc="warpedtemplate",
            space="corobl",
            hemi="{hemi}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rm {params.template_anat} {output.warped}  -r  {input.warp} {params.xfm_corobl}"
