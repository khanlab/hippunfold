import math


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
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.coords} -threshold {params.threshold_tofrom} 1 0 {input.mask} -multiply -o {output}"


rule register_midthickness_syn:
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
                desc="syn",
                **inputs.subj_wildcards,
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                dir="IO",
                label="{label}",
                suffix="invxfm.nii.gz",
                to="{inout}",
                space="corobl",
                hemi="{hemi}",
                desc="syn",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    threads: 16
    conda:
        "../envs/ants.yaml"
    params:
        metric="MeanSquares",
        metric_weight=1,
        radius=0,
        convergence="50x50x50",
        convergence_thresh="1e-6",
        convergence_window=10,
        shrink_factors="4x2x1",
        smoothing_sigmas="2x1x0vox",
        syn_gradient_step=0.1,
        syn_update_field_variance=config["inner_outer_reg_smoothing"],
        syn_total_field_variance=0,
    log:
        bids_log(
            "register_midthickness_ants",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            to="{inout}",
        ),
    shadow:
        "minimal"
    shell:
        r"""
        (
            tmp_prefix=ants_midthickness

            ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} antsRegistration \
                --dimensionality 3 \
                --float 1 \
                --interpolation Linear \
                --use-histogram-matching 0 \
                --transform SyN[{params.syn_gradient_step},{params.syn_update_field_variance},{params.syn_total_field_variance}] \
                --metric {params.metric}[{input.fixed},{input.moving},{params.metric_weight},{params.radius}] \
                --convergence [{params.convergence},{params.convergence_thresh},{params.convergence_window}] \
                --shrink-factors {params.shrink_factors} \
                --smoothing-sigmas {params.smoothing_sigmas} \
                --output [${{tmp_prefix}},${{tmp_prefix}}Warped.nii.gz] \
                --verbose 1

            mv ${{tmp_prefix}}0Warp.nii.gz {output.warp}
            mv ${{tmp_prefix}}0InverseWarp.nii.gz {output.invwarp}
        ) &> {log}
        """


rule register_midthickness_greedy:
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
    params:
        update_field_sigma=math.sqrt(float(config["inner_outer_reg_smoothing"])),
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
                desc="greedy",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    threads: 16
    conda:
        "../envs/greedy.yaml"
    log:
        bids_log(
            "register_midthickness",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            to="{inout}",
        ),
    shell:
        "greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} -n 50x50x50 -s {params.update_field_sigma}vox 0.707vox -o {output.warp} &> {log}"


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
            desc=config["inner_outer_reg_method"],
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
    conda:
        "../envs/greedy.yaml"
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
            desc="{desc}",
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
                    desc="{desc}",
                    **inputs.subj_wildcards,
                )
            )
        ),
    group:
        "subj"
    conda:
        "../envs/workbench.yaml"
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule warp_midthickness_to_inout:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
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
            desc=config["inner_outer_reg_method"],
            **inputs.subj_wildcards,
        ),
    output:
        surf_gii=temp(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname,inner|outer}.surf.gii",
                den="native",
                desc="nostruct",
                space="corobl",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
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
