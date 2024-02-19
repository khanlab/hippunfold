def get_smoothing_opt(wildcards):
    """sets the smoothness of the greedy template shape injection deformation"""

    gradient_sigma = 1.732 * float(config["template_seg_smoothing_factor"])
    warp_sigma = 0.7071 * float(config["template_seg_smoothing_factor"])

    return f"-s {gradient_sigma}vox {warp_sigma}vox"


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
        return wildcards


rule template_reg:
    input:
        fixed_img=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{modality}.nii.gz".format(
                modality=get_modality_suffix(config["modality"])
            ),
            space="corobl",
            desc="preproc",
            hemi="{hemi}"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        moving_img=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]][
            get_modality_suffix(config["modality"])
        ].format(**flipped(wildcards)),
        xfm_corobl=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
        general_opts="-d 3 -m NCC 2x2x2",
        smoothing_opts=get_smoothing_opt,
        iteration_opts="-n 100x50x10",  #default -n 100x100
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
            hemi=hemi_constrained_wildcard,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} "
        " {params.smoothing_opts} {params.iteration_opts} "
        " -i {input.fixed_img} {params.moving_img} -it {params.xfm_corobl} -o {output.warp}"


rule warp_template_dseg:
    input:
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
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        template_dseg=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["dseg"].format(
            **flipped(wildcards)
        ),
        interp_opt="-ri LABEL 0.2vox",
    output:
        inject_seg=bids(
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
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {params.template_dseg} {output.inject_seg}  -r {input.warp}"


rule warp_template_coords:
    input:
        template_dir=Path(download_dir) / "template" / config["template"],
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
        template_coords=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["coords"].format(
            **flipped(wildcards)
        ),
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
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {params.template_coords} {output.init_coords}  -r {input.warp}"


rule warp_template_anat:
    input:
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
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        template_anat=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]][config["modality"]].format(
            **flipped(wildcards)
        ),
        xfm_corobl=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
    output:
        warped=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix=f"{config['modality']}.nii.gz",
            desc="warpedtemplate",
            space="corobl",
            hemi=hemi_constrained_wildcard,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} -rf {input.ref} -rm {params.template_anat} {output.warped}  -r  {input.warp} {params.xfm_corobl}"


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
            hemi="{hemi,L|R}",
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input.nii} -flip x -popas FLIPPED "
        " {input.unflip_ref} -push FLIPPED -copy-transform -o {output.nii} "
