
rule import_manualseg_with_resample:
    input:
        in_img=partial(get_single_bids_input, component="manualseg"),
    params:
        resample_cmd=(
            ""
            if config["resample_manualseg"] == None
            else "-resample {res}".format(res=config["resample_manualseg"])
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
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -int 0 {params.resample_cmd} {params.crop_cmd} -o {output}"
