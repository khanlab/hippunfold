# ruleorder: resample_metric_to_atlas > atlas_metric_to_unfold_nii
# import glob.glob


rule slice3d_to_2d:
    # this is needed so antsMultivariateTemplateConstruction2 will believe the data is truly 2d
    input:
        img=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        img=bids(
            root=work,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfoldeven2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    run:
        import nibabel as nib

        img = nib.load(input.img)
        matrix = img.get_fdata()[:, :, 0]
        nib.Nifti1Image(matrix, img.affine).to_filename(output.img)


rule write_image_pairs_csv:
    input:
        metric_nii=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=work,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfoldeven2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=["gyrification", "curvature", "thickness"],
            hemi=wildcards.hemi,
            label=wildcards.label,
        ),
    output:
        images_csv="template/pairs_{hemi}_{label}.csv",
    group:
        "subj"
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
        images_csv="template/pairs_{hemi}_{label}.csv",
    params:
        num_modalities=len(
            config["atlas_files"][config["gen_template_name"]]["metric_wildcards"]
        ),
        warp_prefix="template/warp_hemi-{hemi}_label-{label}_",
    output:
        warp_example="template/warp_hemi-{hemi}_label-{label}_template0.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["autotop"]  # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        "../envs/ants.yaml"
    shell:
        "antsMultivariateTemplateConstruction2.sh -d 2 -o {params.warp_prefix} -n 0 -l 0 -k {params.num_modalities} {input}"


# def warp_names(wildcards,input):
#     # TODO: this currently won't work with wildcards.session
#     warp_file=glob(f"template/warp_hemi-{hemi}_label-{label}_*_sub-{wildcards.subject}_*1Warp.nii.gz")
#     return warp_file

# # Now we can plug into unfold_reg.smk to out everything in space-unfolreg
# rule mv_unfold_reg:
#     input:
#         warp_example="template/warp_hemi-{hemi}_label-{label}_template0.nii.gz",
#     params:
#         filename = warp_names
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
