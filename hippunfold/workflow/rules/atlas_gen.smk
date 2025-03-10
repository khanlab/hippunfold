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
    group:
        "subj"
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


### final atlas files ###


ruleorder: gen_atlas_surfs > download_extract_atlas
ruleorder: avg_metrics > download_extract_atlas
ruleorder: maxprob_subfields > download_extract_atlas


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
        midthickness_surf=directory(Path(download_dir) / "atlas" / "{atlas}")
        / config["atlas_files"]["{atlas}"]["surf_gii"],
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("pyvista")
    script:
        "../scripts/surf_gen.py"


rule gen_atlas_densities:
    input:
        midthickness_surf=directory(Path(download_dir) / "atlas" / "{atlas}")
        / config["atlas_files"]["{atlas}"]["surf_gii"],
    output:
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
    shell:
        "cp {input} {output}"  # TODO get actual densities


rule avg_metrics:
    input:
        metric_nii=get_metric_template,
        surf=expand(
            directory(Path(download_dir) / "atlas" / config["gen_template_name"])
            + "/"
            + config["atlas_files"]["mytemplate"]["surf_gii"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
        ),
    output:
        metric=expand(
            directory(Path(download_dir) / "atlas" / config["gen_template_name"])
            + "/"
            + config["atlas_files"]["mytemplate"]["metric_gii"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
            metric=config["atlas_files"]["mytemplate"]["metric_wildcards"],
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    shell:
        "c3d {input.metric_nii} -mean tmp.nii.gz && "
        "wb_command -volume-to-surface-mapping tmp.nii.gz {input.surf} {output} -nearest"


"""
rule maxprob_subfields:
    input:
        in_img=inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="subfields.label.gii",
                space="unfold",
                den="{density}",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
            atlas=config["gen_template_name"],
            density="unfoldiso",
        ),
    output:
        label_gii=directory(Path(download_dir) / "atlas" / "{atlas}")
        / config["atlas_files"]["{atlas}"]["label_gii"],
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    shell:
        "c3d {input} -vote -type uchar -o {output}"
"""
