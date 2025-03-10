from glob import glob


rule slice3d_to_2d:
    # this is needed so antsMultivariateTemplateConstruction2 will believe the data is truly 2d
    input:
        img=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        img=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    run:
        import nibabel as nib

        img = nib.load(input.img)
        matrix = img.get_fdata()[:, :, 0]
        nib.Nifti1Image(matrix, img.affine).to_filename(output.img)


def get_cmd_templategen_subj_csv(wildcards, input, output):

    cmds = []

    # write to csv
    csv_line = ",".join(input.metric_nii)
    cmds.append(f"echo {csv_line} > {output.metrics_csv}")

    return " && ".join(cmds)


rule templategen_subj_csv:
    input:
        metric_nii=expand(
            bids(
                root=work,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=["gyrification", "curvature", "thickness"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_templategen_subj_csv,
    output:
        metrics_csv=bids(
            root=work,
            datatype="anat",
            suffix="metrics.csv",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    shell:
        "{params.cmd}"


def get_inputs_templategen_combined_csv(wildcards):

    if "hemi" in inputs[config["modality"]].zip_lists:
        # hemi is an input wildcard,
        #  so it will be already included when we expand
        expand_hemi = {}
    else:
        # hemi is not an input wildcard,
        # so we additionally expand using the config hemi
        expand_hemi = {"hemi": config["hemi"]}

    return inputs[config["modality"]].expand(
        bids(
            root=work,
            datatype="anat",
            suffix="metrics.csv",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        label=wildcards.label,
        **expand_hemi,
    )


rule template_gen_combined_csv:
    input:
        metrics_csvs=get_inputs_templategen_combined_csv,
    params:
        cmd=lambda wildcards, input, output: f"cat {input.metrics_csvs} > {output.metrics_csv}",
    output:
        metrics_csv="template/concat_label-{label}_metrics.csv",
    shell:
        "{params.cmd}"


rule gen_atlas_reg_ants:
    input:
        images_csv="template/concat_label-{label}_metrics.csv",
    params:
        num_modalities=len(
            config["atlas_files"][config["gen_template_name"]]["metric_wildcards"]
        ),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
        multires="-f 6x4 -s 3x2 -q 50x20",  #only two low-res stages for now, to speed up for debugging workflow..
    #        multires=" -f 6x4x2x1 -s 3x2x1x0 -q 100x100x70x20 ",
    output:
        avgtemplate_dir=directory("template/avgtemplate_{label}"),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]  # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        "../envs/ants.yaml"
    shell:
        "antsMultivariateTemplateConstruction2.sh "
        " {params.multires} "
        " -d 2 -o {params.warp_prefix} -n 0 -l 0 -k {params.num_modalities} {input.images_csv} "


rule rename_output_warp:
    input:
        metric_nii=expand(
            bids(
                root=work,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=["gyrification", "curvature", "thickness"],
            allow_missing=True,
        ),
        avgtemplate_dir="template/avgtemplate_{label}",
    params:
        glob_input_warp=lambda wildcards, input: "{ants_prefix}input*-{filename}-1Warp.nii.gz".format(
            ants_prefix=rules.gen_atlas_reg_ants.params.warp_prefix,
            filename=Path(input.metric_nii[0]).name.removesuffix(".nii.gz"),
        ),
        glob_input_invwarp=lambda wildcards, input: "{ants_prefix}input*-{filename}-1InverseWarp.nii.gz".format(
            ants_prefix=rules.gen_atlas_reg_ants.params.warp_prefix,
            filename=Path(input.metric_nii[0]).name.removesuffix(".nii.gz"),
        ),
    output:
        warp=bids(
            root=work,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="unfold",
            to="avgunfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=work,
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="unfold",
            to="avgunfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    shell:
        "cp {params.glob_input_warp} {output.warp} && "
        "cp {params.glob_input_invwarp} {output.invwarp}"


def get_metric_template(wildcards, input):

    num_modalities = (
        len(config["atlas_files"][config["gen_template_name"]]["metric_wildcards"]),
    )

    warp_prefix = "template/warp_label-{label}_"
    matrics = sorted(glob(f"{warp_prefix}_template?.nii.gz"))
    return matrics


rule gen_atlas_surfs:
    input:
        metric_nii=get_metric_template,
    params:
        z_level=get_unfold_z_level,
    output:
        surf=config["atlas_files"]["mytemplate"]["surf_gii"],
    script:
        "../scripts/surf_gen.py"


rule avg_metrics:
    input:
        metric_nii=get_metric_template,
    output:
        metric=config["atlas_files"]["mytemplate"]["metric_gii"],
    shell:
        "c3d {input} -mean tmp.nii.gz && "
        "wb_command -volume-to-surface-mapping tmp.nii.gz {output} -trilinear"


rule maxprob_subfields:
    input:
        in_img=partial(get_single_bids_input, component="dsegsubfields"),
    output:
        maxprob_subfields=config["atlas_files"]["mytemplate"]["label_gii"],
    shell:
        "c3d {input} -vote -type uchar -o {output}"
