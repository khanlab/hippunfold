
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
            metric=config["atlas_metrics"],
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
        metrics_csv=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
        ),
    shell:
        "{params.cmd}"


rule gen_atlas_reg_ants:
    input:
        metrics_csv=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
        ),
    params:
        num_modalities=len(config["atlas_metrics"]),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
        multires="-f 6x4 -s 3x2 -q 50x20",  #only two low-res stages for now, to speed up for debugging workflow..
    #        multires=" -f 6x4x2x1 -s 3x2x1x0 -q 100x100x70x20 ",
    output:
        avgtemplate_dir=directory(
            bids_atlas(
                root=get_atlas_dir(),
                template=config["new_atlas_name"],
                label="{label}",
                suffix="antstemplate",
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]  # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        "../envs/ants.yaml"
    shell:
        "antsMultivariateTemplateConstruction2.sh "
        " {params.multires} "
        " -d 2 -o {params.warp_prefix} -n 0 -l 0 -k {params.num_modalities} {input.metrics_csv} "


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
        avgtemplate_dir=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            suffix="antstemplate",
        ),
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
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=work,
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
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
        avgtemplate_dir=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            suffix="antstemplate",
        ),
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
            i=config["atlas_metrics"].index(wildcards.metric),
        ),
    output:
        metric=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            suffix="{metric}.nii.gz",
        ),
    script:
        "../scripts/set_metric_nii_header.py"


rule gen_atlas_surfs:
    input:
        metric_nii=expand(
            bids_atlas(
                root=get_atlas_dir(),
                template=config["new_atlas_name"],
                label="{label}",
                suffix="{metric}.nii.gz",
            ),
            metric=config["atlas_metrics"],
            allow_missing=True,
        ),
    params:
        z_level=get_unfold_z_level,
    output:
        surf_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            suffix="{surfname,midthickness|inner|outer}.surf.gii",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("pyvista")
    script:
        "../scripts/surf_gen.py"


rule avgtemplate_metric_vol_to_surf:
    input:
        metric_nii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            suffix="{metric}.nii.gz",
        ),
        midthickness=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            suffix="{metric}.shape.gii",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    shell:
        "wb_command -volume-to-surface-mapping {input.metric_nii} {input.midthickness} {output.metric_gii} -trilinear"


# TODO: label voting and mesh resampling..


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
            atlas=config["new_atlas_name"],
            density="unfoldiso",
        ),
    output:
        label_gii=directory(get_atlas_dir() / "{atlas}")
        / config["atlas_files"]["{atlas}"]["label_gii"],
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    shell:
        "c3d {input} -vote -type uchar -o {output}"
"""
