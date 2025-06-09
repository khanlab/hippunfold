import pandas as pd


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
                root=root,
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
        metrics_csv=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="metrics.csv",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    shell:
        "{params.cmd}"


rule template_gen_combined_csv:
    input:
        metrics_csvs=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
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
            root=root,
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
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
        ),
    params:
        num_modalities=len(config["atlas_metrics"]),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
    output:
        avgtemplate_dir=directory(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                suffix="antstemplate",
            )
        ),
        # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        "../envs/ants.yaml"
    shell:
        "antsMultivariateTemplateConstruction2.sh "
        " -d 2 -o {params.warp_prefix} -n 0 -l 0 -k {params.num_modalities} {input.metrics_csv} "


rule copy_avgtemplate_warps:
    input:
        metric_nii=expand(
            bids(
                root=root,
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
            root=root,
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
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                suffix="invwarp.nii.gz",
                from_="unfold",
                to=config["new_atlas_name"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    shell:
        "cp {params.glob_input_warp} {output.warp} && "
        "cp {params.glob_input_invwarp} {output.invwarp}"


rule copy_avgtemplate_metric:
    """ copy avgtemplate metric out of avgtemplate folder"""
    input:
        avgtemplate_dir=bids_atlas(
            root=root,
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
                root=root,
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
                root=root,
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
            root=root,
            template=config["new_atlas_name"],
            desc="badhdr",
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric,[a-zA-Z0-9]+}.nii.gz",
        ),
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/set_metric_nii_header.py"


rule reset_header_2d_warp_atlasgen:
    """ adjusts header to match the original data
     (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfold2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            metric=config["atlas_metrics"][0],
            **expand_hemi(),
        )[0],
        nii=bids(
            root=root,
            datatype="warps",
            suffix="{warpsuffix}.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="warps",
            suffix="{warpsuffix,warp|invwarp}.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            desc="3D",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/set_metric_nii_header.py"


rule make_metric_ref:
    input:
        metric_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz".format(metric=config["atlas_metrics"][0]),
        ),
    output:
        metric_ref=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                resample="{resample}",
                suffix="metricref.nii.gz",
            )
        ),
    conda:
        "../envs/c3d.yaml"
    shell:
        "c2d {input} -resample {wildcards.resample}% -scale 0 -shift 1 -binarize -o {output}"


rule gen_unfold_atlas_mesh:
    input:
        metric_ref=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            resample="{resample}",
            suffix="metricref.nii.gz",
        ),
    params:
        z_level=get_unfold_z_level,
    output:
        surf_gii=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                resample="{resample}",
                space="unfold",
                hemi="{hemi}",
                suffix="{surfname,midthickness|inner|outer}.surf.gii",
            )
        ),
    conda:
        "../envs/pyvista.yaml"
    script:
        "../scripts/gen_unfold_atlas_mesh.py"


def get_unfold_mesh_resample(wildcards):

    ind_resample = config["density_choices"].index(wildcards.density)
    resample = config["resample_factors"][ind_resample]

    return bids_atlas(
        root=root,
        template=config["new_atlas_name"],
        label="{label}",
        resample="{resample}",
        space="unfold",
        hemi="{hemi}",
        suffix="{surfname}.surf.gii",
    ).format(resample=resample, **wildcards)


rule copy_unfold_mesh_resample_to_density:
    input:
        get_unfold_mesh_resample,
    output:
        surf_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            space="unfold",
            hemi="{hemi}",
            suffix="{surfname,midthickness|inner|outer}.surf.gii",
        ),
    shell:
        "cp {input} {output}"


rule avgtemplate_metric_vol_to_surf:
    input:
        metric_nii=bids_atlas(
            root=root,
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
            den="{density}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            hemi="{hemi}",
            suffix="{metric}.shape.gii",
        ),
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -volume-to-surface-mapping {input.metric_nii} {input.midthickness} {output.metric_gii} -trilinear"


rule warp_subj_unfold_surf_to_avg:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            desc="3D",
            **inputs.subj_wildcards,
        ),
    params:
        cmd=get_cmd_warp_surface_2d_warp,
    output:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfoldavg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    shadow:
        "minimal"
    shell:
        "{params.cmd}"


rule resample_subj_native_surf_to_avg:
    input:
        subj_native=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        subj_unfold=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="unfoldavg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        atlas_unfold=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            resample="{resample}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        native_resampled=bids(
            root=root,
            datatype="surf",
            label="{label}",
            space="corobl",
            resample="{resample}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -surface-resample {input.subj_native} {input.subj_unfold} {input.atlas_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check"


rule warp_subfields_to_avg:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        img=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="subfields.nii.gz",
                space=config["new_atlas_name"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/ants.yaml"
    shell:
        "antsApplyTransforms -d 2 -i {input.img} -o {output.img} -t {input.warp} -r {input.img} -n NearestNeighbor -v"


rule vote_subfield_labels:
    input:
        subfield_niis=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
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
                root=root,
                template=config["new_atlas_name"],
                desc="subfields",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
                label="{label}",
            )
        ),
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/majority_voting.py"


rule reset_header_2d_subfields_nii:
    """ adjusts header to match the original data
     (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
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
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfields",
            suffix="dseg.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsfixhdr",
            suffix="dseg.nii.gz",
        ),
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/set_metric_nii_header.py"


rule import_avg_subfields_as_label:
    """read in volumetric subfield dseg, used for group_create_atlas only"""
    input:
        vol_dseg=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsfixhdr",
            suffix="dseg.nii.gz",
        ),
        label_list=Path(workflow.basedir) / "../resources/atlas-v2/labellist_withdg.txt",
    output:
        label_dseg=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldswithlbl",
            suffix="dseg.nii.gz",
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -volume-label-import {input.vol_dseg} {input.label_list} {output.label_dseg}"


rule avgtemplate_subfield_voted_vol_to_surf:
    input:
        subfields_nii=bids_atlas(
            root=root,
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
            den="{density}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            hemi="{hemi}",
            suffix="dseg.label.gii",
        ),
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -volume-label-to-surface-mapping {input.subfields_nii} {input.midthickness} {output.metric_gii}"


# input function for the rule aggregate
def get_atlas_inputs(wildcards):

    files = []

    for label in config["autotop_labels"]:
        for hemi in config["hemi"]:

            files.extend(
                expand(
                    bids_atlas(
                        root=get_atlas_dir(),
                        template=config["new_atlas_name"],
                        label=label,
                        hemi=hemi,
                        den="{density}",
                        suffix="{metric}.shape.gii",
                    ),
                    metric=config["atlas_metrics"],
                    density=config["density_choices"],
                )
            )
            files.extend(
                expand(
                    bids_atlas(
                        root=get_atlas_dir(),
                        template=config["new_atlas_name"],
                        label=label,
                        hemi=hemi,
                        den="{density}",
                        space="unfold",
                        suffix="{surfname}.surf.gii",
                    ),
                    surfname=["inner", "outer", "midthickness"],
                    density=config["density_choices"],
                )
            )
            if label == "hipp":
                files.extend(
                    expand(
                        bids_atlas(
                            root=get_atlas_dir(),
                            template=config["new_atlas_name"],
                            label="hipp",
                            hemi=hemi,
                            den="{density}",
                            suffix="dseg.label.gii",
                        ),
                        density=config["density_choices"],
                    )
                )

    return files


rule write_template_json:
    input:
        get_atlas_inputs,
    params:
        template_description={
            "Identifier": config["new_atlas_name"],
            "metric_wildcards": config["atlas_metrics"],
            "label_wildcards": config["autotop_labels"],
            "hemi_wildcards": config["hemi"],
            "density_wildcards": config["density_choices"],
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
    script:
        "../scripts/write_template_json.py"
