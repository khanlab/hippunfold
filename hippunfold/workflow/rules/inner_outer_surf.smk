

rule compute_halfthick_mask:
    input:
        coords=lambda wildcards: bids(
            root=root,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="coords.nii.gz",
            desc=config["laminar_coords_method"],
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        mask=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        threshold_tofrom=lambda wildcards: (
            "0.5 1" if wildcards.inout == "inner" else "0 0.5"
        ),
    output:
        nii=temp(
            bids(
                root=root,
                datatype="coords",
                dir="IO",
                label="{label}",
                suffix="mask.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    shell:
        "c3d {input.coords} -threshold {params.threshold_tofrom} 1 0 {input.mask} -multiply -o {output}"


rule register_midthickness:
    input:
        fixed=bids(
            root=root,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        moving=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        warp=temp(
            bids(
                root=root,
                datatype="warps",
                dir="IO",
                label="{label}",
                suffix="xfm.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 16
    conda:
        conda_env("greedy")
    log:
        bids_log(
            "register_midthickness",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            to="{inout}",
        ),
    shell:
        "greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} -n 30x0 -o {output.warp} &> {log}"


rule apply_halfsurf_warp_to_img:
    input:
        fixed=bids(
            root=root,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        moving=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            dir="IO",
            label="{label}",
            suffix="xfm.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        warped=temp(
            temp(
                bids(
                    root=root,
                    datatype="coords",
                    dir="IO",
                    label="{label}",
                    suffix="warpedmask.nii.gz",
                    to_="{inout}",
                    space="corobl",
                    hemi="{hemi}",
                    **inputs.subj_wildcards,
                )
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("greedy")
    shell:
        "greedy -d 3  -rf {input.fixed} -rm {input.moving} {output.warped}  -r {input.warp} "


# TODO: rename the warps here according to existing custom (type=itk -> type=ras)
rule convert_inout_warp_from_itk_to_world:
    input:
        warp=bids(
            root=root,
            datatype="warps",
            dir="IO",
            label="{label}",
            suffix="xfm.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        warp=temp(
            temp(
                bids(
                    root=root,
                    datatype="warps",
                    dir="IO",
                    label="{label}",
                    suffix="xfmras.nii.gz",
                    to="{inout}",
                    space="corobl",
                    hemi="{hemi}",
                    **inputs.subj_wildcards,
                )
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule warp_midthickness_to_inout:
    input:
        surf_gii=bids(
            root=root,
            datatype="surfnative",
            suffix="midthickness.surf.gii",
            den="native",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            dir="IO",
            label="{label}",
            suffix="xfmras.nii.gz",
            to="{surfname}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    output:
        surf_gii=temp(
            bids(
                root=root,
                datatype="surfnative",
                suffix="{surfname,inner|outer}.surf.gii",
                den="native",
                space="corobl",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("workbench")
    shadow:
        "minimal"
    group:
        "subj"
    log:
        bids_log(
            "warp_midthickness_to_inout",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            to="{surfname}",
        ),
    shell:
        """
        (
            wb_command -volume-to-surface-mapping {input.warp} {input.surf_gii} warp.shape.gii -trilinear &&
            wb_command -surface-coordinates-to-metric {input.surf_gii} coords.shape.gii &&
            wb_command -metric-math 'COORDS + WARP' warpedcoords.shape.gii -var COORDS coords.shape.gii -var WARP warp.shape.gii &&
            wb_command -surface-set-coordinates {input.surf_gii} warpedcoords.shape.gii {output.surf_gii}
        ) &> {log}
        """
