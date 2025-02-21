ruleorder: resample_metric_to_atlas > atlas_metric_to_unfold_nii


rule write_image_pairs_csv:
    input:
        metric_nii=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        images_csv = expand(
            "template/pairs_{hemi}_{label}.csv",
            hemi=config["hemi"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
        ),
    run:
        with open(output.images_csv, "w") as f:
            subjects = sorted(set(re.search(r"sub-(\w+)", i).group(1) for i in input.metric_nii))
            for subj in subjects:
                metrics = [i for i in input.metric_nii if f"sub-{subj}" in i]
                f.write(",".join(metrics) + "\n")

rule gen_atlas_reg_ants:
    input:
        images_csv = expand(
            "template/pairs_{hemi}_{label}.csv",
            hemi=config["hemi"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
        ),
    params:
        num_modalities = len(congif["atlas_files"][config["gen_template_name"]]["metric_wildcards"])
    output:
        warp = bids(
            root=template,
            preffix="warp",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
        ),
    shell:
        "antsMultivariateTemplateConstruction2.sh -d 2 -o {output} -n 0 -l 0 -k {params.num_modalities} {input}"