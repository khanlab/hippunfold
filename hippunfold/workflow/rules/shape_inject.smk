# Template-based segmentation supports templates that have only a single hemisphere
# by flipping it


def get_input_splitseg_for_shape_inject(wildcards):
    if config["modality"] == "dsegtissue":
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dsegsplit",
            space="corobl",
            hemi="{hemi}",
        ).format(**wildcards)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dsegsplit",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
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

    conda:
        conda_env("c3d")
    shell:
        "mkdir -p {output} && "
        "c3d {input} -retain-labels {params.labels} -split -foreach -smooth {params.smoothing_stdev} -endfor -oo {output}/label_%02d.nii.gz"


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


rule resample_template_dseg_tissue_for_reg:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    params:
        resample_cmd="-resample-mm {res}".format(
            res=config["resample_dseg_for_templatereg"]
        ),
        crop_cmd="-trim 5vox",  #leave 5 voxel padding
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissueresampled",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -int 0 {params.resample_cmd} {params.crop_cmd} -o {output}"


rule template_shape_reg:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissueresampled",
            hemi="{hemi}",
            suffix="dsegsplit",
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
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
    group:
        "subj"

    conda:
        conda_env("greedy")
    threads: 8
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            suffix="templateshapereg.txt",
        ),
    shell:
        #affine (with moments), then greedy
        "greedy -threads {threads} {params.general_opts} {params.affine_opts} {params.img_pairs} -o {output.matrix}  &> {log} && "
        "greedy -threads {threads} {params.general_opts} {params.greedy_opts} {params.img_pairs} -it {output.matrix} -o {output.warp} &>> {log}"


rule dilate_dentate_pd_src_sink:
    """ The PD src/sink labels can disappear after label propagation
    as they are very small. This dilates them into relative background labels"""
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    params:
        src_label=config["laplace_labels"]["dentate"]["PD"]["src"][0],
        sink_label=config["laplace_labels"]["dentate"]["PD"]["sink"][0],
        src_bg=2,
        sink_bg=10,
        struc_elem_size=3,
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissuedilated",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    group:
        "subj"
    conda:
        conda_env("neurovis")
    script:
        "../scripts/dilate_dentate_pd_src_sink.py"


rule template_shape_inject:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="{label}tissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            label="{label}",
            space="corobl",
            hemi="{hemi}",
        ),
        matrix=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        interp_opt="-ri LABEL 0.1mm",  # smoothing sigma = 100micron
    output:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="inject",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
        ),
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            suffix="templateshapeinject.txt",
            hemi="{hemi}",
            label="{label}",
        ),
    group:
        "subj"

    conda:
        conda_env("greedy")
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.upsampled_ref} -rm {input.template_seg} {output.inject_seg}  -r {input.warp} {input.matrix} &> {log}"


rule inject_init_laplace_coords:
    """ TODO: this may not be needed anymore """
    input:
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
        ),
        matrix=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
        coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="init",
            space="template",
            hemi="{hemi}",
        ),
    params:
        interp_opt="-ri LIN",
    output:
        init_coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi}",
        ),
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="{label}",
            suffix="injectcoords.txt",
            desc="init",
            hemi="{hemi}",
        ),
    group:
        "subj"

    conda:
        conda_env("greedy")
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.upsampled_ref} -rm {input.coords} {output.init_coords}  -r {input.warp} {input.matrix} &> {log}"


rule reinsert_subject_labels:
    """ c3d command to:
                1) get the labels to retain
                2) reslice to injected seg
                3) set injected seg to zero where retained labels are
                4) and add this to retained labels"""
    input:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="inject",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
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
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
        ),
    group:
        "subj"

    conda:
        conda_env("c3d")
    shell:
        "c3d {input.subject_seg} -retain-labels {params.labels} -popas LBL "
        " -int 0 {input.inject_seg} -as SEG -push LBL -reslice-identity -popas LBL_RESLICE "
        "-push LBL_RESLICE -threshold 0 0 1 0 -push SEG -multiply "
        "-push LBL_RESLICE -add -o {output.postproc_seg}"
