def get_smoothing_opt(wildcards):
    """sets the smoothness of the greedy template shape injection deformation"""

    gradient_sigma = 1.732 * float(config["template_seg_smoothing_factor"])
    warp_sigma = 0.7071 * float(config["template_seg_smoothing_factor"])

    return f"-s {gradient_sigma}vox {warp_sigma}vox"


rule template_reg:
    input:
        fixed_img=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        moving_img=bids(
            root=root,
            datatype="anat",
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            desc="template",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        general_opts="-d 3 -m NCC 2x2x2",
        smoothing_opts=get_smoothing_opt,
        iteration_opts="-n 100x50x10",  #default -n 100x100
    output:
        warp=temp(
            bids(
                root=root,
                **inputs.subj_wildcards,
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="greedytemplatereg",
                from_="template",
                to="subject",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/greedy.yaml"
    log:
        bids_log("template_reg", **inputs.subj_wildcards, hemi="{hemi}"),
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} "
        " {params.smoothing_opts} {params.iteration_opts} "
        " -i {input.fixed_img} {input.moving_img} -o {output.warp} &> {log}"


rule warp_template_dseg:
    input:
        upsampled_ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=root,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
        template_dseg=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="hipptissue",
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    params:
        interp_opt="-ri LABEL 0.2vox",
    output:
        inject_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="postproc",
                space="corobl",
                hemi="{hemi}",
                label="hipp",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/greedy.yaml"
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.upsampled_ref} -rm {input.template_dseg} {output.inject_seg}  -r {input.warp}"


rule warp_template_dseg_dentate:
    input:
        upsampled_ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=root,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
        template_dseg=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="dentatetissue",
            space="template",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    params:
        interp_opt="-ri LABEL 0.2vox",
    output:
        inject_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="postproc",
                space="corobl",
                hemi="{hemi}",
                label="dentate",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/greedy.yaml"
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.upsampled_ref} -rm {input.template_dseg} {output.inject_seg}  -r {input.warp}"
