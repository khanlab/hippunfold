ruleorder: resample_metric_to_atlas > atlas_metric_to_unfold_nii


rule write_image_pairs_csv:
    input:
        metric_nii=expand(bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        )),
    output:
        images_csv=expand(
            "template/pairs_{hemi}_{label}.csv",
            hemi="{hemi}",
            label="{label}",
        ),
    run:
        with open(output.images_csv, "w") as f:
            subjects = sorted(
                set(re.search(r"sub-(\w+)", i).group(1) for i in input.metric_nii)
            )
            for subj in subjects:
                metrics = [i for i in input.metric_nii if f"sub-{subj}" in i]
                f.write(",".join(metrics) + "\n")


rule gen_atlas_reg_ants:
    input:
        images_csv=expand(
            "template/pairs_{hemi}_{label}.csv",
            hemi=config["hemi"],
            label=config["atlas_files"]["mytemplate"]["label_wildcards"],
        ),
    params:
        num_modalities=len(
            congif["atlas_files"][config["gen_template_name"]]["metric_wildcards"]
        ),
    output:
        warp=expand(
            "template/warp_hemi-{hemi}_label-{label}_"
            hemi="{hemi}",
            label="{label}",
        ),
    shell:
        "antsMultivariateTemplateConstruction2.sh -d 2 -o {output} -n 0 -l 0 -k {params.num_modalities} {input}"

def warp_names(wildcards,input):
    # TODO: get the actual filename output by ants template creation
    warp_file=f"template/warp_hemi-{hemi}_label-{label}_sub-{}"
    return warp_file

# # Now we can plug into unfold_reg.smk to out everything in space-unfolreg
# rule mv_unfold_reg:
#     input:
#         warp=bids(
#             root=template,
#             preffix="warp",
#             space="unfold",
#             hemi="{hemi}",
#             label="{label}",
#         ),
#     params:
#         filename = ""
#     output:
#         warp=bids(
#             root=work,
#             suffix="xfm.nii.gz",
#             datatype="warps",
#             desc="{desc}",
#             from_="{from}",
#             to="{to}",
#             space="{space}",
#             type_="itk2d",
#             hemi="{hemi}",
#             label="{label}",
#             **inputs.subj_wildcards,
#         ),
#     shell:
#         "mv {params.filename} {output.warp}"

# rule: gen_atlas_surfs
# # this should output to our atlas directory. Should also generate variious densities (decimation to target?)

# rule: avg_metrics:
# # this should output to our atlas directory

# rule: maxprob_subfields
# # this should output to our atlas directory

