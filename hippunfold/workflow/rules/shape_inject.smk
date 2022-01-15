def get_input_for_shape_inject(wildcards):
    if config["modality"] == "cropseg":
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}"
        ).format(**wildcards)
    elif get_modality_key(config["modality"]) == "seg":
        modality_suffix = get_modality_suffix(config["modality"])
        seg = (
            bids(
                root=work,
                datatype="anat",
                **config["subj_wildcards"],
                suffix="dseg.nii.gz",
                space="corobl",
                hemi="{hemi}",
                from_="{modality_suffix}"
            ).format(**wildcards, modality_suffix=modality_suffix),
        )
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}"
        ).format(**wildcards)
    return seg


def get_input_splitseg_for_shape_inject(wildcards):
    if config["modality"] == "cropseg":
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dsegsplit",
            space="corobl",
            hemi="{hemi}"
        ).format(**wildcards)

    elif get_modality_key(config["modality"]) == "seg":
        modality_suffix = get_modality_suffix(config["modality"])
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dsegsplit",
            space="corobl",
            hemi="{hemi}",
            from_="{modality_suffix}"
        ).format(**wildcards, modality_suffix=modality_suffix)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dsegsplit",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}"
        ).format(**wildcards)
    return seg


rule prep_segs_for_greedy:
    input:
        "{prefix}_dseg.nii.gz",
    params:
        labels=" ".join(str(label) for label in config["shape_inject"]["labels_reg"]),
        smoothing_stdev=config["shape_inject"]["label_smoothing_stdev"],
    output:
        directory("{prefix}_dsegsplit"),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "mkdir -p {output} && "
        "c3d {input} -retain-labels {params.labels} -split -foreach -smooth {params.smoothing_stdev} -endfor -oo {output}/label_%02d.nii.gz"


rule import_template_shape:
    input:
        template_seg=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "tpl-upenn",
            "tpl-upenn_desc-hipptissue_dseg.nii.gz",
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **config["subj_wildcards"],
            desc="hipptissue",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


def get_image_pairs(wildcards, input):
    """This rule requires snakemake 6.4.0, since it uses the new feature to execute if input files are not found"""

    import errno

    if not os.path.exists(input.subject_seg):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), input.subject_seg
        )
    if not os.path.exists(input.template_seg):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), input.template_seg
        )

    args = []
    # prep_segs_for_greedy creates label_{i} images for each entry in labels_reg,
    # but the numbering will be from 1 to N (not the numbers in the list)
    for label in range(1, len(config["shape_inject"]["labels_reg"]) + 1):
        subject_label = f"{input.subject_seg}/label_{label:02d}.nii.gz"
        template_label = f"{input.template_seg}/label_{label:02d}.nii.gz"

        if not os.path.exists(subject_label):
            print(f"Warning: {subject_label} does not exist, not using in registration")
            continue
        if not os.path.exists(template_label):
            print(
                f"Warning: {template_label} does not exist, not using in registration"
            )
            continue

        args.append("-i")
        args.append(subject_label)  # subject is fixed
        args.append(template_label)  # template is moving
    return " ".join(args)


def get_inject_scaling_opt(wildcards):
    """sets the smoothness of the greedy template shape injection deformation"""

    gradient_sigma = 1.732 * float(config["inject_template_smoothing_factor"])
    warp_sigma = 0.7071 * float(config["inject_template_smoothing_factor"])

    return f"-s {gradient_sigma}vox {warp_sigma}vox"


rule template_shape_reg:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **config["subj_wildcards"],
            desc="hipptissue",
            suffix="dsegsplit"
        ),
        subject_seg=get_input_splitseg_for_shape_inject,
    params:
        general_opts="-d 3 -m SSD",
        affine_opts="-moments 2 -det 1",
        greedy_opts=get_inject_scaling_opt,
        img_pairs=get_image_pairs,
    output:
        matrix=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi,Lflip|R}"
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
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
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            hemi="{hemi,Lflip|R}",
            suffix="templateshapereg.txt"
        ),
    shell:
        #affine (with moments), then greedy
        "greedy -threads {threads} {params.general_opts} {params.affine_opts} {params.img_pairs} -o {output.matrix}  &> {log} && "
        "greedy -threads {threads} {params.general_opts} {params.greedy_opts} {params.img_pairs} -it {output.matrix} -o {output.warp} &>> {log}"


rule template_shape_inject:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **config["subj_wildcards"],
            desc="hipptissue",
            suffix="dseg.nii.gz"
        ),
        subject_seg=get_input_for_shape_inject,
        matrix=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}"
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
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
            desc="inject",
            space="corobl",
            hemi="{hemi,Lflip|R}"
        ),
    log:
        bids(
            root='logs',
            **config["subj_wildcards"],
            suffix="templateshapeinject.txt",
            hemi="{hemi,Lflip|R}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.subject_seg} -rm {input.template_seg} {output.inject_seg}  -r {input.warp} {input.matrix} &> {log}" 


rule inject_init_laplace_coords:
    input:
        coords=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "tpl-upenn",
            "tpl-upenn_dir-{dir}_label-{autotop}_coords.nii.gz",
        ),
        subject_seg=get_input_for_shape_inject,
        matrix=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}"
        ),
        warp=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
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
            hemi="{hemi}"
        ),
    log:
        bids(
            root='logs',
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="injectcoords.txt",
            desc="init",
            hemi="{hemi}"
        ),        
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.subject_seg} -rm {input.coords} {output.init_coords}  -r {input.warp} {input.matrix} &> {log}"


rule reinsert_subject_labels:
    input:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="inject",
            space="corobl",
            hemi="{hemi}"
        ),
        subject_seg=get_input_for_shape_inject,
    params:
        labels=" ".join(
            str(label) for label in config["shape_inject"]["labels_reinsert"]
        ),
    output:
        postproc_seg=bids(
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
    shell:
        "c3d {input.subject_seg} -retain-labels {params.labels} -popas LBL -push LBL -threshold 0 0 1 0 {input.inject_seg} -multiply -push LBL -add -o {output.postproc_seg}"


rule unflip_postproc:
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
        "c3d {input} -flip x {output}"
