# get _b500.nii.gz from the bids dataset - this is a preprocessed avg b500 image


rule resample_hippdwi_to_template:
    """ Hipp DWI is already corobl, but just needs to be cropped
    to get L and R subvolumes. We use predefined X and Y 
    bounding boxes, and keep all Z slices (since is Z 
    is smaller than in template). """
    input:
        b500=inputs["hippb500"].path,
    params:
        resample_dim=config["hippdwi_opts"]["resample_dim"],
        bbox_x=lambda wildcards: config["hippdwi_opts"]["bbox_x"][wildcards.hemi],
        bbox_y=config["hippdwi_opts"]["bbox_y"],
    output:
        crop_b500=bids(
            root=work,
            datatype="dwi",
            hemi="{hemi,L|R}",
            space="corobl",
            suffix="b500.nii.gz",
            **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -resample {params.resample_dim} -as UPSAMPLED "
        " -push UPSAMPLED -cmv -pop -popas COORDY -popas COORDX "
        " -push COORDX -thresh {params.bbox_x} 1 0 -as MASKX "
        " -push COORDY -thresh {params.bbox_y} 1 0 -as MASKY "
        " -push MASKX -push MASKY -multiply "
        " -push UPSAMPLED -multiply "
        " -trim 0vox -o {output}"


rule cp_b500_to_anat_dir:
    input:
        nii=bids(
            root=work,
            datatype="dwi",
            **inputs.subj_wildcards,
            suffix="b500.nii.gz",
            space="corobl",
            hemi="{hemi}"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            desc="preproc",
            suffix="hippb500.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"
