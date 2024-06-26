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
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
        ),
    output:
        nii=bids(
            root=work,
            datatype="coords",
            space="unfold",
            label="{autotop}",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
        ),
        unfold_phys_coords_nii=bids(
            root=work,
            space="unfold",
            label="hipp",
            datatype="coords",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
        ),
        coords_io=get_laminar_coords,
        native_ref_coords_nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="hipp",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"]
        ),
        unfold_phys_coords_nii=bids(
            root=work,
            space="unfold",
            label="dentate",
            datatype="coords",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
        ),
        native_ref_coords_nii=bids(
            root=work,
            datatype="coords",
            space="corobl",
            hemi="{hemi}",
            label="dentate",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
            hemi="{hemi}",
            suffix="create_warps-dentate.txt"
        ),
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/create_warps.py"


rule expand_unfolded_warps:
    """unfolded space registration in 2D expanded to 3D"""
    input:
        warp2d=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN",
            from_="{from}",
            to="{to}",
            space="unfold",
            type_="itk",
            hemi="{hemi}"
        ),
        unfold_phys_coords_nii=bids(
            root=work,
            space="unfold",
            label="hipp",
            datatype="coords",
            suffix="refcoords.nii.gz",
            **config["subj_wildcards"]
        ),
    output:
        warp3d=bids(
            root=work,
            **config["subj_wildcards"],
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="SyN3D",
            from_="{from}",
            to="{to}",
            space="unfold",
            type_="itk",
            hemi="{hemi}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/expand_2Dwarp.py"


def get_unfold2unfoldatlas(wildcards):
    if config["no_unfolded_reg"]:
        fn = []
    else:
        fn = (
            bids(
                root=work,
                **config["subj_wildcards"],
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="SyN3D",
                from_="subject",
                to=config["atlas"][0],
                space="unfold",
                type_="itk",
                hemi="{hemi}"
            ),
        )
    return fn


rule compose_warps_native_to_unfold:
    """ Compose warps from native to unfold """
    input:
        unfold2unfoldatlas=get_unfold2unfoldatlas,
        corobl2unfold=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="corobl",
            to="unfold",
            mode="image"
        ),
        ref=bids(
            root=root,
            space="unfold",
            label="{autotop}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **config["subj_wildcards"]
        ),
        native2corobl=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
    output:
        bids(
            root=root,
            datatype="warps",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="{native_modality}",
            to="unfold",
            mode="image"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="composexfm.txt",
            hemi="{hemi}",
            from_="{native_modality}",
            to="unfold",
            mode="image"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "ComposeMultiTransform 3 {output} -R {input.ref} {input.unfold2unfoldatlas} {input.corobl2unfold} {input.native2corobl} &> {log}"


def get_unfoldatlas2unfold(wildcards):
    if config["no_unfolded_reg"]:
        fn = []
    else:
        fn = (
            bids(
                root=work,
                **config["subj_wildcards"],
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="SyN3D",
                from_=config["atlas"][0],
                to="subject",
                space="unfold",
                type_="itk",
                hemi="{hemi}",
            ),
        )
    return fn


def get_cmd_compose_warps_unfold_to_crop_native(wildcards, input, output):
    if config["no_unfolded_reg"]:
        cmd = f"antsApplyTransforms -o [{output.unfold2cropnative},1] -r {input.ref} -t [{input.native2corobl},1] -t {input.unfold2corobl} -i {input.unfold_ref} -v"
    else:
        cmd = f"antsApplyTransforms -o [{output.unfold2cropnative},1] -r {input.ref} -t [{input.native2corobl},1] -t {input.unfold2corobl} -t  {input.unfoldatlas2unfold} -i {input.unfold_ref} -v"
    return cmd


rule compose_warps_unfold_to_crop_native:
    """ Compose warps from unfold to crop native """
    input:
        unfoldatlas2unfold=get_unfoldatlas2unfold,
        unfold2corobl=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="corobl",
            mode="image"
        ),
        ref=bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{native_modality}",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
        native2corobl=bids(
            root=work,
            datatype="warps",
            **config["subj_wildcards"],
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="itk"
        ),
        unfold_ref=bids(
            root=root,
            space="unfold",
            label="{autotop}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **config["subj_wildcards"]
        ),
    output:
        unfold2cropnative=bids(
            root=root,
            datatype="warps",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="xfm.nii.gz",
            hemi="{hemi}",
            from_="unfold",
            to="{native_modality}",
            mode="image"
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            label="{autotop}",
            suffix="composexfm.txt",
            hemi="{hemi}",
            from_="unfold",
            to="{native_modality}",
            mode="image"
        ),
    params:
        cmd=get_cmd_compose_warps_unfold_to_crop_native,
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "{params.cmd}  &> {log}"
