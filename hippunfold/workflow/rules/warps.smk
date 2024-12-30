def reg_to_template_cmd(wildcards, input, output):
    ref = str(
        Path(input.template_dir)
        / config["template_files"][config["template"]][wildcards.modality].format(
            **wildcards
        ),
    )
    if config["no_reg_template"]:
        cmd = f"reg_resample -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {input.xfm_identity}; cp {input.xfm_identity} {output.xfm_ras}"
    elif config["rigid_reg_template"]:
        cmd = f"reg_aladin -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {output.xfm_ras} -rigOnly"
    else:
        cmd = f"reg_aladin -flo {input.flo} -ref {ref} -res {output.warped_subj} -aff {output.xfm_ras}"
    return cmd


rule reg_to_template:
    """ generic for T1w or T2w right now """
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="preproc",
            suffix="{modality}.nii.gz"
        ),
        xfm_identity=os.path.join(workflow.basedir, "..", config["xfm_identity"]),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        cmd=reg_to_template_cmd,
    output:
        warped_subj=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="{modality,T1w|T2w}.nii.gz",
            space=config["template"],
            desc="affine"
        ),
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality,T1w|T2w}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            suffix="reg.txt",
            from_="{modality,T1w|T2w}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}"


rule convert_template_xfm_ras2itk:
    input:
        bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{reg_suffix}",
            to=config["template"],
            desc="affine",
            type_="ras"
        ),
    output:
        bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{reg_suffix}",
            to=config["template"],
            desc="affine",
            type_="itk"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input}  -oitk {output}"


# now have subject -> template transform, can compose that with template -> corobl to get subject -> corobl
rule compose_template_xfm_corobl:
    input:
        sub_to_std=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to=config["template"],
            desc="affine",
            type_="itk"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        std_to_cor=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["xfm_corobl"].format(
            **wildcards
        ),
    output:
        sub_to_cor=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input.sub_to_std} -itk {params.std_to_cor} -mult -oitk {output}"


rule invert_template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    output:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="T1w",
            to="corobl",
            desc="affineInverse",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -inv -o {output}"


rule template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    output:
        xfm_ras=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality,T1w|T2w}",
            to="corobl",
            desc="affine",
            type_="ras"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -o {output}"


rule create_native_coord_ref:
    input:
        coords_ap=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    output:
        nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -cmp -omc {output}"


# unfold ref nifti
rule create_unfold_ref:
    params:
        dims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["dims"]
        ),
        voxdims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["voxdims"]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.autotop]["origin"]
        ),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.autotop]["orient"],
    output:
        nii=bids(
            root=root,
            space="unfold",
            label="{autotop}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}"


# this was unfold_phys_coords.nii in matlab implementation
rule create_unfold_coord_map:
    input:
        nii=bids(
            root=root,
            space="unfold",
            label="{autotop}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
    output:
        nii=bids(
            root=work,
            datatype="coords",
            space="unfold",
            label="{autotop}",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.nii} -cmp -omc {output.nii}"


def get_laminar_coords(wildcards):
    if "laplace" in config["laminar_coords_method"]:
        coords_io = bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        )
    elif "equivolume" in config["laminar_coords_method"]:
        coords_io = bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="hipp",
            suffix="coords.nii.gz",
            desc="equivol",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        )
    return coords_io


rule create_warps_hipp:
    input:
        unfold_ref_nii=bids(
            root=root,
            space="unfold",
            label="hipp",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
        unfold_phys_coords_nii=bids(
            root=work,
            space="unfold",
            label="hipp",
            datatype="coords",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
        coords_ap=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        coords_pd=bids(
            root=work,
            datatype="coords",
            dir="PD",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        coords_io=get_laminar_coords,
        native_ref_coords_nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: config["laplace_labels"]["AP"]["gm"],
    resources:
        mem_mb=16000,
    output:
        warp_unfold2native=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="hipp",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="surface"
        ),
        warp_native2unfold=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="hipp",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="surface"
        ),
        warpitk_unfold2native=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="hipp",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="image"
        ),
        warpitk_native2unfold=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="hipp",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="image"
        ),
    group:
        "subj"
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            suffix="create_warps-hipp.txt"
        ),
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/create_warps.py"


rule create_warps_dentate:
    input:
        unfold_ref_nii=bids(
            root=root,
            space="unfold",
            label="dentate",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards
        ),
        unfold_phys_coords_nii=bids(
            root=work,
            space="unfold",
            label="dentate",
            datatype="coords",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
        coords_ap=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="dentate",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        coords_pd=bids(
            root=work,
            datatype="coords",
            dir="PD",
            label="dentate",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        coords_io=bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="dentate",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        native_ref_coords_nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="dentate",
            suffix="refcoords.nii.gz",
            **inputs.subj_wildcards
        ),
        labelmap=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: config["laplace_labels"]["PD"]["sink"],
    resources:
        mem_mb=16000,
    output:
        warp_unfold2native=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="dentate",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="surface"
        ),
        warp_native2unfold=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="dentate",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="surface"
        ),
        warpitk_unfold2native=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="dentate",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="image"
        ),
        warpitk_native2unfold=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            label="dentate",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="image"
        ),
    group:
        "subj"
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            suffix="create_warps-dentate.txt"
        ),
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/create_warps.py"
