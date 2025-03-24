
rule slice_3d_to_2d:
    """This is needed so antsMultivariateTemplateConstruction2 will believe the data is truly 2d"""
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
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("neurovis")
    script:
        "../scripts/slice_3d_to_2d.py"


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


rule template_gen_combined_csv:
    input:
        metrics_csvs=lambda wildcards: inputs[config["modality"]].expand(
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
        ),
    params:
        cmd=lambda wildcards, input, output: f"cat {input.metrics_csvs} > {output.metrics_csv}",
    output:
        metrics_csv=bids_atlas(
            root=work,
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
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
        ),
    params:
        num_modalities=len(config["atlas_metrics"]),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
        # multires="-f 6x4 -s 3x2 -q 50x20",  #only two low-res stages for now, to speed up for debugging workflow..
        multires=" -f 6x4x2x1 -s 3x2x1x0 -q 100x100x70x20 ",
    output:
        avgtemplate_dir=directory(
            bids_atlas(
                root=work,
                template=config["new_atlas_name"],
                label="{label}",
                suffix="antstemplate",
            )
        ),
    container:
        config["singularity"]["autotop"]  # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        conda_env("ants")
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
            metric=config["atlas_metrics"],
            allow_missing=True,
        ),
        avgtemplate_dir=bids_atlas(
            root=work,
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
    shell:
        "cp {params.glob_input_warp} {output.warp} && "
        "cp {params.glob_input_invwarp} {output.invwarp}"


rule copy_avgtemplate_metric:
    """ copy avgtemplate metric out of avgtemplate folder"""
    input:
        avgtemplate_dir=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            suffix="antstemplate",
        ),
    params:
        in_metric=lambda wildcards, input: "{avgtemplate_dir}/template{i}.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            i=config["atlas_metrics"].index(wildcards.metric),
        ),
    output:
        metric=temp(
            bids_atlas(
                root=work,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                desc="badhdr",
                suffix="{metric}.nii.gz",
            )
        ),
    shell:
        "cp {params.in_metric} {output.metric}"


rule reset_header_2d_metric_nii:
    """ adjusts header to match the original data
     (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=lambda wildcards: inputs[config["modality"]].expand(
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
        nii=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            desc="badhdr",
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric,[a-zA-Z0-9]+}.nii.gz",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("neurovis")
    script:
        "../scripts/set_metric_nii_header.py"


rule gen_atlas_surfs:
    input:
        metric_nii=expand(
            bids_atlas(
                root=work,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
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
            hemi="{hemi}",
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
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz",
        ),
        midthickness=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.shape.gii",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    shell:
        "wb_command -volume-to-surface-mapping {input.metric_nii} {input.midthickness} {output.metric_gii} -trilinear"


rule warp_subfields_to_avg:
    input:
        img=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
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
    output:
        img=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("ants")
    shell:
        "antsApplyTransforms -d 2 -i {input.img} -o {output.img} -t {input.warp} -r {input.img} -n NearestNeighbor -v"


rule vote_subfield_labels:
    input:
        subfield_niis=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=work,
                datatype="anat",
                suffix="subfields.nii.gz",
                space=config["new_atlas_name"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            **expand_hemi(),
        ),
    output:
        subfields_voted=temp(
            bids_atlas(
                root=work,
                template=config["new_atlas_name"],
                desc="subfields",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
                label="{label}",
            )
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("neurovis")
    script:
        "../scripts/majority_voting.py"


rule reset_header_2d_subfields_nii:
    """ adjusts header to match the original data
     (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=lambda wildcards: inputs[config["modality"]].expand(
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
            metric=config["atlas_metrics"],
            **expand_hemi(),
        )[0],
        nii=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfields",
            suffix="dseg.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsfixhdr",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("neurovis")
    script:
        "../scripts/set_metric_nii_header.py"


rule import_avg_subfields_as_label:
    """read in volumetric subfield dseg, used for group_create_atlas only"""
    input:
        vol_dseg=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsfixhdr",
            suffix="dseg.nii.gz",
        ),
        label_list=Path(workflow.basedir) / "../resources/atlas-v2/labellist_withdg.txt",
    output:
        label_dseg=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldswithlbl",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    group:
        "subj"
    shell:
        "wb_command -volume-label-import {input.vol_dseg} {input.label_list} {output.label_dseg}"


rule avgtemplate_subfield_voted_vol_to_surf:
    input:
        subfields_nii=bids_atlas(
            root=work,
            template=config["new_atlas_name"],
            desc="subfieldswithlbl",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
            label="{label}",
        ),
        midthickness=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="dseg.label.gii",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    shell:
        "wb_command -volume-label-to-surface-mapping {input.subfields_nii} {input.midthickness} {output.metric_gii}"


rule write_template_json:
    input:
        metrics=expand(
            bids_atlas(
                root=get_atlas_dir(),
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                suffix="{metric}.shape.gii",
            ),
            hemi=config["hemi"],
            label=config["autotop_labels"],
            metric=config["atlas_metrics"],
        ),
        surfs=expand(
            bids_atlas(
                root=get_atlas_dir(),
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                space="unfold",
                suffix="{surfname}.surf.gii",
            ),
            hemi=config["hemi"],
            label=config["autotop_labels"],
            surfname=["inner", "outer", "midthickness"],
        ),
        labels=expand(
            bids_atlas(
                root=get_atlas_dir(),
                template=config["new_atlas_name"],
                label="hipp",
                hemi="{hemi}",
                suffix="dseg.label.gii",
            ),
            hemi=config["hemi"],
        ),
    params:
        template_description={
            "Identifier": config["new_atlas_name"],
            "metric_wildcards": config["atlas_metrics"],
            "label_wildcards": config["autotop_labels"],
            "hemi_wildcards": config["hemi"],
            "Authors": [""],
            "Acknowledgements": "",
            "BIDSVersion": "",
            "HowToAcknowledge": "",
            "License": "MIT",
            "Name": "HippUnfold surface-based template and subfield atlas",
            "RRID": "",
            "ReferencesAndLinks": [""],
            "TemplateFlowVersion": "",
        },
    output:
        json=str(
            Path(
                bids(
                    root=get_atlas_dir(),
                    tpl=config["new_atlas_name"],
                )
            )
            / "template_description.json"
        ),
