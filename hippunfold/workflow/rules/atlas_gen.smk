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


#for template generation, need line for each subj/hemi
#  will do this by making a file for each subject/hemi first, then concatenating those

def get_cmd_templategen_subj_csv(wildcards, input, output):
 
    cmds=[]
    
#    #copy nifti to new filenames (for ants naming)
#    for in_nii,out_nii in zip(input.metric_nii,output.metric_nii):
#        cmds.append(f"cp {in_nii} {out_nii}")

    #write to csv
    csv_line=",".join(input.metric_nii)
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
            metric=["gyrification", "curvature", "thickness"], #TODO: use config
            allow_missing=True,
        )
    params:
        cmd=get_cmd_templategen_subj_csv
    output:
        metrics_csv=
            bids(
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
            **expand_hemi)



rule template_gen_combined_csv:
    input:
        metrics_csvs=get_inputs_templategen_combined_csv
    params:
        cmd=lambda wildcards, input, output:
            f"cat {input.metrics_csvs} > {output.metrics_csv}"
    output:
        metrics_csv="template/concat_label-{label}_metrics.csv", #TODO: use bids()
    shell:
        "{params.cmd}"


def get_cmd_rename_output_warps(wildcards, input):
 
    num_modalities=len(
            config["atlas_files"][config["gen_template_name"]]["metric_wildcards"]
        ),
   
    warp_prefix="template/warp_label-{label}_"

    warps=sorted(glob(f'{warp_prefix}input*-1Warp.nii.gz'))
    invwarps=sorted(glob(f'{warp_prefix}input*-1InverseWarp.nii.gz'))
    
    cmds=[]
    
    for i,(warp,invwarp) in enumerate(zip(warps,invwarps)):
        cmds.append(f"mv {warp} {warp_prefix}subindex-{i}_warp.nii.gz")
        cmds.append(f"mv {invwarp} {warp_prefix}subindex-{i}_invwarp.nii.gz")

    return " && ".join(cmds)



rule gen_atlas_reg_ants:
    input:
        images_csv="template/concat_label-{label}_metrics.csv", #TODO: use bids()
    params:
        num_modalities=len(
            config["atlas_files"][config["gen_template_name"]]["metric_wildcards"]
        ),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
    output:
        avgtemplate_dir=directory("template/avgtemplate_{label}")
    group:
        "subj"
    container:
        config["singularity"]["autotop"]  # note antsMultivariateTemplateConstruction2 is not in the container right now!
    conda:
        "../envs/ants.yaml"
    shell:
        "antsMultivariateTemplateConstruction2.sh "
        " -f 6x4 -s 3x2 -q 50x20 " #only two low-res stages for now, to speed up for debugging workflow..
        #" -f 6x4x2x1 -s 3x2x1x0 -q 100x100x70x20 " 
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
            metric=["gyrification", "curvature", "thickness"][0], #TODO: use config
            allow_missing=True,
        ),
        avgtemplate_dir="template/avgtemplate_{label}"
    params:
        glob_input_warp=lambda wildcards, input: f"{input.avgtemplate_dir}/input*-{Path(input.metric_nii[0]).name.removesuffix('.nii.gz')}-1Warp.nii.gz",
        glob_input_invwarp=lambda wildcards, input: f"{input.avgtemplate_dir}/input*-{Path(input.metric_nii[0]).name.removesuffix('.nii.gz')}-1InverseWarp.nii.gz",
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



rule warp_surf_to_avg_template:
    """ added this to verify the warped surfaces --
    NOTE: this isn't working currently, since the warps are now strictly 2D """
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
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
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="avgunfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shadow:
        "minimal"
    shell:
        "wb_command -volume-to-surface-mapping {input.warp} {input.surf_gii} warp.shape.gii -trilinear && "
        "wb_command -surface-coordinates-to-metric {input.surf_gii} coords.shape.gii && "
        "wb_command -metric-math 'COORDS + WARP' warpedcoords.shape.gii -var COORDS coords.shape.gii -var WARP warp.shape.gii && "
        "wb_command -surface-set-coordinates  {input.surf_gii} warpedcoords.shape.gii {output.surf_gii}"






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
