
rule import_dseg_tissue:
    input:
        in_img=partial(get_single_bids_input, component="dsegtissue"),
    params:
        resample_cmd=(
            ""
            if config["resample_dsegtissue"] == None
            else "-resample {res}".format(res=config["resample_dsegtissue"])
        ),
        crop_cmd="-trim 5vox",  #leave 5 voxel padding
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi,L|R}",
        ),

    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -int 0 {params.resample_cmd} {params.crop_cmd} -o {output}"
