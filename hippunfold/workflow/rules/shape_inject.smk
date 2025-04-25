# Template-based segmentation supports templates that have only a single hemisphere
# by mapping it to the flipped version of the other hemisphere.
# If a template has both L and R files, then we set hemi_constrained_wildcard to L|R.
# If a hemisphere is missing data, then we set it to flip that, e.g. if L missing, then use Lflip|R

hemi_constraints = []
if config["template"] in config["template_based_segmentation"]:
    for hemi in config["hemi"]:
        if hemi in config["template_based_segmentation"][config["template"]]["hemi"]:
            hemi_constraints.append(hemi)
        else:
            hemi_constraints.append(f"{hemi}flip")

hemi_constrained_wildcard = "{{hemi,{constraints}}}".format(
    constraints="|".join(hemi_constraints)
)


def flipped(wildcards):
    """function to map hemi in wildcards from Lflip to R, or Rflip to L,
    for use in rules where e.g. the output wildcard is Lflip, but for the input, R is desired, such as
    when mapping a R hemi dseg to the Lflip hemisphere of a subject."""

    if wildcards.hemi == "L" or wildcards.hemi == "R":
        return wildcards
    elif wildcards.hemi == "Lflip":
        wildcards.hemi = "R"
        return wildcards
    elif wildcards.hemi == "Rflip":
        wildcards.hemi = "L"
        return wildcard


def get_input_for_shape_inject(wildcards):
    if config["modality"] == "cropseg":
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            space="corobl",
            hemi="{hemi}",
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
                from_="{modality_suffix}",
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
            hemi="{hemi}",
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
            hemi="{hemi}",
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
            from_="{modality_suffix}",
        ).format(**wildcards, modality_suffix=modality_suffix)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
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
    container:
        config["singularity"]["autotop"]
    shell:
        "mkdir -p {output} && "
        "c3d {input} -retain-labels {params.labels} -split -foreach -smooth {params.smoothing_stdev} -endfor -oo {output}/label_%02d.nii.gz"


rule import_template_shape:
    input:
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        template_seg=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["inject_template"]]["dseg"].format(
            **flipped(wildcards)
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **config["subj_wildcards"],
            desc="hipptissue",
            hemi=hemi_constrained_wildcard,
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {params.template_seg} {output.template_seg}"


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
            hemi="{hemi}",
            suffix="dsegsplit"
        ),
        subject_seg=get_input_splitseg_for_shape_inject,
    params:
        general_opts="-d 3 -m SSD",
        affine_opts="-moments 2",
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
            hemi=hemi_constrained_wildcard,
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
            hemi=hemi_constrained_wildcard,
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
            hemi=hemi_constrained_wildcard,
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
            hemi="{hemi}",
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
            hemi="{hemi}",
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
            hemi="{hemi}",
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
            hemi=hemi_constrained_wildcard,
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            suffix="templateshapeinject.txt",
            hemi=hemi_constrained_wildcard,
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
            hemi="{hemi}",
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
            hemi="{hemi}",
        ),
        template_dir=Path(download_dir) / "template" / config["inject_template"],
    params:
        coords=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["inject_template"]]["coords"].format(
            **wildcards
        ),
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
            hemi=hemi_constrained_wildcard,
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="injectcoords.txt",
            desc="init",
            hemi=hemi_constrained_wildcard,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.subject_seg} -rm {params.coords} {output.init_coords}  -r {input.warp} {input.matrix} &> {log}"


rule unflip_init_coords:
    """Unflip the Lflip init coords"""
    input:
        nnunet_seg=bids(
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
        unflip_ref=get_input_for_shape_inject,
    output:
        nnunet_seg=bids(
            root=work,
            datatype="coords",
            **config["subj_wildcards"],
            dir="{dir}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi,L|R}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nnunet_seg} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nnunet_seg} "


rule reinsert_subject_labels:
    input:
        inject_seg=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="inject",
            space="corobl",
            hemi="{hemi}",
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
            hemi=hemi_constrained_wildcard,
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
        unflip_ref=get_input_for_shape_inject,
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi,L|R}",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nii} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nii} "
