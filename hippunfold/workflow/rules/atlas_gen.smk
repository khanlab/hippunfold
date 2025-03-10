from glob import glob


def expand_hemi():
    if "hemi" in inputs[config["modality"]].zip_lists:
        # hemi is an input wildcard,
        #  so it will be already included when we expand
        return {}
    else:
        # hemi is not an input wildcard,
        # so we additionally expand using the config hemi
        return {"hemi": config["hemi"]}


# TODO: get these from config eventually
reg_metrics = ["gyrification", "curvature", "thickness"]


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
            metric=reg_metrics,
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
        **expand_hemi(),
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
        num_modalities=len(reg_metrics),
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


rule copy_avgtemplate_warps:
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
        glob_input_warp=lambda wildcards, input: "{avgtemplate_dir}/input*-{filename}-1Warp.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            filename=Path(input.metric_nii[0]).name.removesuffix(".nii.gz"),
        ),
        glob_input_invwarp=lambda wildcards, input: "{avgtemplate_dir}/input*-{filename}-1InverseWarp.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
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


rule copy_avgtemplate_metric:
    """ copy avgtemplate metric, adjusting header to match the original data
     (since this seems to get garbled in z by ants)"""
    input:
        avgtemplate_dir="template/avgtemplate_{label}",
        ref_metric=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=work,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            metric=wildcards.metric,
            **expand_hemi(),
        )[0],
    params:
        in_metric=lambda wildcards, input: "{avgtemplate_dir}/template{i}.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            i=reg_metrics.index(wildcards.metric),
        ),
    output:
        metric="template/avgtemplate_label-{label}_{metric}.nii.gz",
    script:
        "../scripts/set_metric_nii_header.py"


def get_metric_template(wildcards, input):

    num_modalities = len(reg_metrics)

    warp_prefix = "template/warp_label-{label}_"
    metrics = sorted(glob(f"{warp_prefix}_template?.nii.gz"))
    return metrics


rule gen_atlas_surfs:
    input:
        metric_nii=expand(
            "template/avgtemplate_label-{label}_{metric}.nii.gz",
            metric=reg_metrics,
            allow_missing=True,
        ),
    params:
        z_level=get_unfold_z_level,
    output:
        surf_gii="template/avgtemplate_space-unfold_label-{label}_{surfname,midthickness|inner|outer}.surf.gii",
    script:
        "../scripts/surf_gen.py"


rule avgtemplate_metric_vol_to_surf:
    input:
        metric_nii="template/avgtemplate_label-{label}_{metric}.nii.gz",
        midthickness="template/avgtemplate_space-unfold_label-{label}_midthickness.surf.gii",
    output:
        metric_gii="template/avgtemplate_label-{label}_{metric}.shape.gii",
    shell:
        "wb_command -volume-to-surface-mapping {input.metric_nii} {input.midthickness} {output.metric_gii} -trilinear"


"""
rule maxprob_subfields:
    input:
        in_img=partial(get_single_bids_input, component="dsegsubfields"),
    output:
        maxprob_subfields=config["atlas_files"]["mytemplate"]["label_gii"],
    shell:
        "c3d {input} -vote -type uchar -o {output}"
"""
