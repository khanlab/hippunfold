rule template_reg:
    input:
        moving_img=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]][config["modality"]],
        ),
        xfm_corobl=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]]["xfm_corobl"],
        ),
        fixed_img=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
    params:
        general_opts="-d 3 -m NCC 2x2x2",
        greedy_opts="-n 100x50x10 -s 4vox 2vox",  #default is 1.732vox, 0.7071vox
    output:
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi,Lflip|R}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} {params.greedy_opts} "
        " -i {input.fixed_img} {input.moving_img} -it {input.xfm_corobl} -o {output.warp}"


rule warp_template_dseg:
    input:
        template_seg=os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]]["dseg"],
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
    params:
        interp_opt="-ri LABEL 0.2vox",
    output:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi,Lflip|R}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {input.template_seg} {output.inject_seg}  -r {input.warp}"


rule warp_template_coords:
    input:
        template_coords=os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]]["coords"],
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
    params:
        interp_opt="-ri NN",
    output:
        init_coords=bids(
            root=work,
            datatype="coords",
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi,R|Lflip}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {input.template_coords} {output.init_coords}  -r {input.warp}"


rule warp_template_anat:
    input:
        template_anat=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]][config["modality"]],
        ),
        ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
        xfm_corobl=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["template_files"][config["template"]]["xfm_corobl"],
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedytemplatereg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}"
        ),
    output:
        warped=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            desc="warpedtemplate",
            space="corobl",
            hemi="{hemi,Lflip|R}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rm {input.template_anat} {output.warped}  -r  {input.warp} {input.xfm_corobl}"


rule unflip_template_dseg:
    input:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}flip",
            **config["subj_wildcards"]
        ),
        unflip_ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi,L}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nii} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nii} "


rule unflip_template_coords:
    input:
        nii=bids(
            root=work,
            datatype="coords",
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi}flip"
        ),
        unflip_ref=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            space="corobl",
            desc="preproc",
            hemi="{hemi}",
        ),
    output:
        nii=bids(
            root=work,
            datatype="coords",
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi,L}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nii} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nii} "
