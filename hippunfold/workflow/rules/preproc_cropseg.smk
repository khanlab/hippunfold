rule import_cropseg:
    input:
        config["input_path"]["cropseg"],
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi,L|R}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


rule lr_flip_seg:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}"
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi,L}flip"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -flip x -o  {output}"
