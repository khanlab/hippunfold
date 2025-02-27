# This is one big rules file for now, but could split it up according to the
# dividers ( --- ) that describe each conceptual part of the pipeline

# --- parameters - will put these in config eventually

surf_thresholds = {"inner": 0, "outer": 1, "midthickness": 0.5}


ruleorder: resample_native_surf_to_std_density > cp_template_to_unfold
ruleorder: atlas_label_to_unfold_nii > atlas_metric_to_unfold_nii
ruleorder: gen_inner_outer_unfoldeven > warp_native_mesh_to_unfold


# --- isosurface generation ---


rule gen_native_mesh:
    input:
        coords=lambda wildcards: bids(
            root=work,
            datatype="coords",
            dir="IO",
            label="{label}",
            suffix="coords.nii.gz",
            desc=config["laminar_coords_method"],
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        nan_mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="nan",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        sink_mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="sink",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        src_mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="IO",
            desc="src",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        threshold=lambda wildcards: surf_thresholds[wildcards.surfname],
        decimate_opts={
            "reduction": 0.05,
            "feature_angle": 25,
            "preserve_topology": True,
        },
        hole_fill_radius=1.0,
        morph_openclose_dist=2,  # mm
        coords_epsilon=0.1,
    output:
        surf_gii=temp(
            bids(
                root=work,
                datatype="surf",
                suffix="{surfname,midthickness}.surf.gii",
                space="corobl",
                desc="nostruct",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    script:
        "../scripts/gen_isosurface.py"


rule update_native_mesh_structure:
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="{space}",
            desc="nostruct",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "cp {input} {output} && wb_command -set-structure {output.surf_gii} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule smooth_surface:
    """ slight smoothing of surface to improve curvature estimation """
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        smoothing_strength=0.8,
        smoothing_iterations=10,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            desc="smoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-smoothing {input} {params.smoothing_strength} {params.smoothing_iterations} {output}"


# --- creating unfold surface from native anatomical, including post-processing


rule get_boundary_vertices:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        label_gii=bids(
            root=work,
            datatype="coords",
            suffix="boundary.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    script:
        "../scripts/get_boundary_vertices.py"


rule map_src_sink_sdt_to_surf:
    """ Maps the distance to src/sink mask """
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        sdt=bids(
            root=work,
            datatype="coords",
            suffix="sdt.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="{srcsink}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        sdt=bids(
            root=work,
            datatype="coords",
            suffix="sdt.shape.gii",
            space="corobl",
            hemi="{hemi}",
            dir="{dir}",
            desc="{srcsink}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -volume-to-surface-mapping {input.sdt} {input.surf_gii} {output.sdt} -trilinear"


rule postproc_boundary_vertices:
    """ ensures non-overlapping and full labelling of AP/PD edges """
    input:
        ap_src=bids(
            root=work,
            datatype="coords",
            suffix="sdt.shape.gii",
            space="corobl",
            hemi="{hemi}",
            dir="AP",
            desc="src",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ap_sink=bids(
            root=work,
            datatype="coords",
            suffix="sdt.shape.gii",
            space="corobl",
            hemi="{hemi}",
            dir="AP",
            desc="sink",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        pd_src=bids(
            root=work,
            datatype="coords",
            suffix="sdt.shape.gii",
            space="corobl",
            hemi="{hemi}",
            dir="PD",
            desc="src",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        pd_sink=bids(
            root=work,
            datatype="coords",
            suffix="sdt.shape.gii",
            space="corobl",
            hemi="{hemi}",
            dir="PD",
            desc="sink",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        edges=bids(
            root=work,
            datatype="coords",
            suffix="boundary.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        min_terminal_vertices=5,  # min number of vertices per src/sink
    output:
        ap=bids(
            root=work,
            datatype="coords",
            suffix="mask.label.gii",
            space="corobl",
            hemi="{hemi}",
            dir="AP",
            desc="srcsink",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        pd=bids(
            root=work,
            datatype="coords",
            suffix="mask.label.gii",
            space="corobl",
            hemi="{hemi}",
            dir="PD",
            desc="srcsink",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    group:
        "subj"
    script:
        "../scripts/postproc_boundary_vertices.py"


rule laplace_beltrami:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        src_sink_mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.label.gii",
            space="corobl",
            hemi="{hemi}",
            dir="{dir}",
            desc="srcsink",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        coords=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    script:
        "../scripts/laplace_beltrami.py"


def get_unfold_z_level(wildcards):
    extent = float(config["unfold_vol_ref"][wildcards.label]["extent"][-1])
    return surf_thresholds[wildcards.surfname] * extent


rule warp_native_mesh_to_unfold:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        coords_AP=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="{label}",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        coords_PD=bids(
            root=work,
            datatype="coords",
            dir="PD",
            label="{label}",
            suffix="coords.shape.gii",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    params:
        vertspace=lambda wildcards: config["unfold_vol_ref"][wildcards.label],
        z_level=get_unfold_z_level,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            desc="nostruct",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyunfold.yaml"
    group:
        "subj"
    script:
        "../scripts/rewrite_vertices_to_flat.py"


rule space_unfold_vertices_evenly:
    """ this irons out the surface to result in more even
        vertex spacing. the resulting shape will be more 
        individual (e.g. the surface area in unfolded space 
        would be similar to native) -- TODO: maybe this is a good 
        way to determine smoothing strenghth and iterations, e.g. 
        use the surface area and vertex spacing as objectives.."""
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            desc="nostruct",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        native_gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        step_size=0.1,
        iterations=1000,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            desc="nostruct",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/pyvista.yaml"
    group:
        "subj"
    log:
        surf_gii=bids(
            root="logs",
            suffix=".txt",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    script:
        "../scripts/space_unfold_vertices.py"


rule gen_inner_outer_unfoldeven:
    input:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            desc="nostruct",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    params:
        z_level=get_unfold_z_level,
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname,inner|outer}.surf.gii",
            desc="nostruct",
            space="unfoldeven",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    run:
        import nibabel as nib

        surf = nib.load(input.surf_gii)
        vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
        vertices[:, 2] = params.z_level
        nib.save(surf, output.surf_gii)


# --- creating inner/outer surfaces from native anatomical (using 3d label deformable registration)


rule compute_halfthick_mask:
    input:
        coords=lambda wildcards: bids(
            root=work,
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
            root=work,
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
                root=work,
                datatype="surf",
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
        "../envs/c3d.yaml"
    shell:
        "c3d {input.coords} -threshold {params.threshold_tofrom} 1 0 {input.mask} -multiply -o {output}"


rule register_midthickness:
    input:
        fixed=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        moving=bids(
            root=work,
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
                root=work,
                datatype="surf",
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
        "../envs/greedy.yaml"
    shell:
        "greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} -n 30x0 -o {output.warp}"


rule apply_halfsurf_warp_to_img:
    input:
        fixed=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="mask.nii.gz",
            to="{inout}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        moving=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            desc="GM",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=work,
            datatype="surf",
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
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="warpedmask.nii.gz",
                to_="{inout}",
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
        "../envs/greedy.yaml"
    shell:
        "greedy -d 3  -rf {input.fixed} -rm {input.moving} {output.warped}  -r {input.warp} "


# TODO: rename the warps here according to existing custom (type=itk -> type=ras)
rule convert_inout_warp_from_itk_to_world:
    input:
        warp=bids(
            root=work,
            datatype="surf",
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
            bids(
                root=work,
                datatype="surf",
                dir="IO",
                label="{label}",
                suffix="xfmras.nii.gz",
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
        "../envs/workbench.yaml"
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule warp_midthickness_to_inout:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=work,
            datatype="surf",
            dir="IO",
            label="{label}",
            suffix="xfmras.nii.gz",
            to="{surfname}",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="ANATOMICAL",
    output:
        surf_gii=bids(
            root=work,
            datatype="surf",
            suffix="{surfname,inner|outer}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        "wb_command -volume-to-surface-mapping {input.warp} {input.surf_gii} warp.shape.gii -trilinear && "
        "wb_command -surface-coordinates-to-metric {input.surf_gii} coords.shape.gii && "
        "wb_command -metric-math 'COORDS + WARP' warpedcoords.shape.gii -var COORDS coords.shape.gii -var WARP warp.shape.gii && "
        "wb_command -surface-set-coordinates  {input.surf_gii} warpedcoords.shape.gii {output.surf_gii}"


# --- affine transforming anatomical surfaces from corobl to other (T1w, T2w) spaces


# warp native surface from corobl to T1w/T2w
rule affine_gii_corobl_to_modality:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=work,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{native_modality}",
            to="corobl",
            desc="affine",
            type_="ras",
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="{native_modality,T1w|T2w}",
            hemi="{hemi}",
            label="{label,hipp|dentate}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}"


# --- calculating metrics on native anatomical surfaces


rule calculate_surface_area:
    input:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-vertex-areas {input} {output}"


rule calculate_legacy_gyrification:
    """new gyrification is ratio of nativearea to unfoldarea (e.g. surface scaling or distortion factor.
    this should be proportional by a constant, to the earlier gyrification on 32k surfaces."""
    input:
        native_surfarea=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        unfold_surfarea=bids(
            root=work,
            datatype="surf",
            suffix="surfarea.shape.gii",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="gyrification.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
        " -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea}"


rule calculate_curvature:
    input:
        gii=bids(
            root=work,
            datatype="surf",
            suffix="midthickness.surf.gii",
            space="corobl",
            desc="smoothed",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="curvature.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-curvature {input} -mean {output}"


rule calculate_thickness:
    input:
        inner=bids(
            root=root,
            datatype="surf",
            suffix="inner.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        outer=bids(
            root=root,
            datatype="surf",
            suffix="outer.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        gii=bids(
            root=root,
            datatype="surf",
            suffix="thickness.shape.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"


# --- resampling using the unfoldreg surface to (legacy) standard densities (0p5mm, 1mm, 2mm, unfoldiso)


def get_unfold_ref_name(wildcards):
    if (
        wildcards.label in config["atlas_files"][config["atlas"]]["label_wildcards"]
        and config["no_unfolded_reg"] == False
    ):
        return "unfoldreg"
    else:
        return "unfold"


def get_unfold_ref(wildcards):
    """function to return either unfoldreg or unfold ref mesh, depending on whether
    unfoldreg can be performed (based on atlas wildcards)"""

    return bids(
        root=root,
        datatype="surf",
        suffix="midthickness.surf.gii",
        space=get_unfold_ref_name(wildcards),
        hemi="{hemi}",
        label="{label}",
        **inputs.subj_wildcards,
    )


rule resample_atlas_subfields_to_std_density:
    """ resamples subfields from a custom atlas mesh (e.g. from an atlas anatomical/native mesh
    warped to unfold space) to a standard density"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        atlas_unfold=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-resample {params.label_gii} {params.atlas_unfold} {input.ref_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


rule resample_native_surf_to_std_density:
    input:
        native=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        native_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="{space,unfoldreg|corobl}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-resample {input.native} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check"


rule resample_native_metric_to_std_density:
    input:
        native_metric=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.{metrictype}.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_unfold=os.path.join(
            workflow.basedir,
            "..",
            "resources",
            "unfold_template_{label}",
            "tpl-avg_space-unfold_den-{density}_midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        metric_resampled=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.{metrictype,shape|func}.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -metric-resample {input.native_metric} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.metric_resampled} -bypass-sphere-check"


rule cp_surf_to_root:
    input:
        native_resampled=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        native_resampled=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="{space}",
            den="{density}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


# --- resampling from atlasnative to native vertices
rule resample_atlas_subfields_to_native_surf:
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        native_unfold=get_unfold_ref,
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        ref_unfold=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="surf",
            suffix="subfields.label.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label,hipp}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-resample {params.label_gii} {params.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


# --- rule for converting atlasnative subfields to unfold nii (e.g. analogous to unfoldiso standard space)
rule atlas_label_to_unfold_nii:
    """converts metric .gii files to .nii - this is a band-aid fix since we make use of the volumetric
        subfield labels in unfold space for painting the native volumetric dseg. Ie this uses the unfold volumetric (or effectively unfoldiso)
        to perform mapping from atlas to subject. better approach would be to adjust the downstream  volumetric dseg function to 
        make use of gifti labels instead.. 
"""
    input:
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
        ref_nii=bids(
            root=work,
            space="unfold",
            label="{label}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards,
        ),
    params:
        label_gii=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["label_gii"].format(**wildcards),
        inner_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="inner", **wildcards
        ),
        outer_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="outer", **wildcards
        ),
        midthickness_surf=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["surf_gii"].format(
            surf_type="midthickness", **wildcards
        ),
    output:
        label_nii=bids(
            root=work,
            datatype="anat",
            suffix="subfields.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {params.label_gii} {params.midthickness_surf} {input.ref_nii} {output.label_nii} "
        " -nearest-vertex 1000"


# --- cifti for native surfs (could use generic rule if we include den wildcard)


def get_inputs_cifti_metric_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_metric"] = (
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="L",
                label="{label}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_metric"] = (
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.shape.gii",
                space="{space}",
                hemi="R",
                label="{label}",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dscalar_metric_cifti_native:
    input:
        unpack(get_inputs_cifti_metric_native),
    params:
        cmd=get_cmd_cifti_metric,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            suffix="{metric}.dscalar.nii",
            space="{space}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmd}"


def get_inputs_cifti_label_native(wildcards):
    files = dict()
    if "L" in config["hemi"]:
        files["left_label"] = (
            bids(
                root=root,
                datatype="surf",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="L",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    if "R" in config["hemi"]:
        files["right_label"] = (
            bids(
                root=root,
                datatype="surf",
                atlas="{atlas}",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="R",
                label="hipp",
                **inputs.subj_wildcards,
            ).format(**wildcards),
        )
    return files


rule create_dlabel_cifti_subfields_native:
    input:
        unpack(get_inputs_cifti_label_native),
    params:
        cmd=get_cmd_cifti_label,
    output:
        cifti=bids(
            root=root,
            datatype="surf",
            atlas="{atlas}",
            suffix="subfields.dlabel.nii",
            space="{space}",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmd}"


# --- create spec file (adapted from rule in gifti.smk)


rule create_spec_file_hipp_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        subfields=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="subfields.label.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        cifti_labels=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="subfields.dlabel.nii",
                atlas="{atlas}",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            atlas=config["atlas"],
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,hipp}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmds}"


rule create_spec_file_dentate_native:
    input:
        metrics=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{metric}.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=get_gifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
        surfs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            surfname=["midthickness"],
            space=["{space}", get_unfold_ref_name(wildcards)],
            allow_missing=True,
        ),
        cifti_metrics=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                suffix="{cifti}.nii",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            cifti=get_cifti_metric_types(wildcards.label),
            allow_missing=True,
        ),
    params:
        cmds=get_cmd_spec_file,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            suffix="surfaces.spec",
            hemi="{hemi,L|R}",
            space="{space}",
            label="{label,dentate}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmds}"


rule merge_lr_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                suffix="surfaces.spec",
                hemi="{hemi}",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            space="{space}",
            suffix="surfaces.spec",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmd}"


rule merge_hipp_dentate_spec_file:
    input:
        spec_files=expand(
            bids(
                root=root,
                datatype="surf",
                suffix="surfaces.spec",
                space="{space}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=config["autotop_labels"],
            allow_missing=True,
        ),
    params:
        cmd=get_cmd_merge_spec,
    output:
        spec_file=bids(
            root=root,
            datatype="surf",
            space="{space}",
            suffix="surfaces.spec",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "{params.cmd}"


rule cp_native_surf_to_root:
    input:
        native=bids(
            root=work,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        native=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"
