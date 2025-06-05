def get_unfold_ref_name(wildcards):
    if (
        wildcards.label in config["atlas_metadata"][config["atlas"]]["label_wildcards"]
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
        den="native",
        hemi="{hemi}",
        label="{label}",
        **inputs.subj_wildcards,
    )


rule cp_atlas_unfold:
    input:
        ref_unfold=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="{surf_name}.surf.gii",
        ),
    output:
        unfold=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name,midthickness}.surf.gii",
            space="unfold",
            den="{density,[0-9k]+}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule resample_native_surf_to_atlas_density:
    input:
        native=bids(
            root=root,
            datatype="surf",
            suffix="{surf_name}.surf.gii",
            space="{space}",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_unfold=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="{surf_name}.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        native_resampled=temp(
            bids(
                root=root,
                datatype="surf",
                suffix="{surf_name,midthickness|inner|outer}.surf.gii",
                space="{space,corobl}",
                den="{density,[0-9k]+}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    log:
        bids_log(
            "resample_native_surf_to_atlas_density",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            space="{space}",
            den="{density}",
            desc="{surf_name}",
        ),
    shell:
        "wb_command -surface-resample {input.native} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check &> {log}"


rule resample_native_metric_to_atlas_density:
    input:
        native_metric=bids(
            root=root,
            datatype="metric",
            suffix="{metric}.{metrictype}.gii",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_unfold=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        metric_resampled=bids(
            root=root,
            datatype="metric",
            suffix="{metric}.{metrictype,shape|func}.gii",
            den="{density,[0-9k]+}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    log:
        bids_log(
            "resample_native_metric_to_atlas_density",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            den="{density,[0-9k]+}",
            desc="{metric}-{metrictype}",
        ),
    shell:
        "wb_command -metric-resample {input.native_metric} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.metric_resampled} -bypass-sphere-check &> {log}"


rule resample_native_coords_to_atlas_density:
    input:
        native_metric=bids(
            root=root,
            datatype="metric",
            dir="{dir}",
            suffix="coords.shape.gii",
            desc="{desc}",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        ref_unfold=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="midthickness.surf.gii",
        ),
        native_unfold=get_unfold_ref,
    output:
        metric_resampled=bids(
            root=root,
            datatype="metric",
            dir="{dir}",
            suffix="coords.shape.gii",
            desc="{desc}",
            den="{density,[0-9k]+}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    log:
        bids_log(
            "resample_native_coords_to_atlas_density",
            **inputs.subj_wildcards,
            hemi="{hemi}",
            label="{label}",
            den="{density,[0-9k]+}",
            dir="{dir}",
            desc="{desc}",
        ),
    shell:
        "wb_command -metric-resample {input.native_metric} {input.native_unfold} {input.ref_unfold} BARYCENTRIC {output.metric_resampled} -bypass-sphere-check &> {log}"


# --- resampling from avgatlas to native vertices
rule resample_atlas_subfields_to_native_surf:
    input:
        native_unfold=get_unfold_ref,
        label_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            den=config["unfoldreg_density"],
            label="{label}",
            suffix="dseg.label.gii",
        ),
        ref_unfold=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            space="unfold",
            den=config["unfoldreg_density"],
            suffix="midthickness.surf.gii",
        ),
    output:
        label_gii=temp(
            bids(
                root=root,
                datatype="metric",
                suffix="subfields.label.gii",
                den="native",
                hemi="{hemi}",
                label="{label,hipp}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-resample {input.label_gii} {input.ref_unfold} {input.native_unfold} BARYCENTRIC {output.label_gii} -bypass-sphere-check"


rule cp_atlas_subfields_label_gii:
    input:
        label_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            den="{density}",
            label="{label}",
            suffix="dseg.label.gii",
        ),
    output:
        label_gii=bids(
            root=root,
            datatype="metric",
            suffix="subfields.label.gii",
            hemi="{hemi}",
            label="{label,hipp}",
            den="{density,[0-9k]+}",
            atlas="{atlas}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule atlas_label_to_unfold_nii:
    """converts metric .gii files to .nii - this is a band-aid fix since we make use of the volumetric
        subfield labels in unfold space for painting the native volumetric dseg. Ie this uses the unfold volumetric (or effectively unfoldiso)
        to perform mapping from atlas to subject. better approach would be to adjust the downstream  volumetric dseg function to 
        make use of gifti labels instead.. 
"""
    input:
        ref_nii=bids(
            root=root,
            space="unfold",
            label="{label}",
            datatype="warps",
            suffix="refvol.nii.gz",
            **inputs.subj_wildcards,
        ),
        label_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            den="{density}",
            label="{label}",
            suffix="dseg.label.gii",
        ),
        midthickness_surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["atlas"],
            hemi="{hemi}",
            label="{label}",
            den="{density}",
            space="unfold",
            suffix="midthickness.surf.gii",
        ),
    output:
        label_nii=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="subfields.nii.gz",
                space="unfold",
                hemi="{hemi}",
                label="{label}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -label-to-volume-mapping {input.label_gii} {input.midthickness_surf} {input.ref_nii} {output.label_nii} "
        " -nearest-vertex 1000"


rule affine_gii_corobl_to_orig:
    input:
        gii=bids(
            root=root,
            datatype="{surfdir}",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            from_="{modality}",
            to="corobl",
            desc="affine",
            type_="ras",
        ),
    output:
        gii=bids(
            root=root,
            datatype="{surfdir}",
            den="{density}",
            suffix="{surfname}.surf.gii",
            space="{modality,T1w|T2w}",
            hemi="{hemi}",
            label="{label,hipp|dentate}",
            **inputs.subj_wildcards,
        ),
    conda:
        "../envs/workbench.yaml"
    group:
        "subj"
    shell:
        "wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}"
