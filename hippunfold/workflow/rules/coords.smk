def get_labels_for_laplace(wildcards):
    if (
        config["skip_inject_template_labels"]
        or config["analysis_level"] == "group_create_atlas"
    ):
        seg = get_input_for_shape_inject(wildcards)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
            label="{label}",
        ).format(**wildcards)
    return seg


def get_gm_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label]["IO"]["gm"]]
    )
    return lbl_list


def get_src_sink_labels(wildcards):
    lbl_list = " ".join(
        [
            str(lbl)
            for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                wildcards.srcsink
            ]
        ]
    )
    return lbl_list


def get_nan_labels(wildcards):
    lbl_list = " ".join(
        [
            str(lbl)
            for lbl in config["laplace_labels"][wildcards.label]["AP"]["sink"]
            + config["laplace_labels"][wildcards.label]["AP"]["src"]
            + config["laplace_labels"][wildcards.label]["PD"]["sink"]
            + config["laplace_labels"][wildcards.label]["PD"]["src"]
        ]
    )
    return lbl_list


rule get_label_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=get_gm_labels,
    output:
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
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


def get_inputs_laplace(wildcards):
    files = dict()
    files["lbl"] = get_labels_for_laplace(wildcards)
    if not config["skip_inject_template_labels"]:
        files["init_coords"] = (
            bids(
                root=work,
                datatype="coords",
                **inputs.subj_wildcards,
                dir="{dir}",
                label="hipp",
                suffix="coords.nii.gz",
                desc="init",
                space="corobl",
                hemi="{hemi}",
            ),
        )
    return files


rule get_src_sink_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=get_src_sink_labels,
    output:
        mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="{srcsink,src|sink}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_src_sink_sdt:
    """calculate signed distance transform (negative inside, positive outside)"""
    input:
        mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
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
            suffix="sdt.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="{srcsink}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -sdt -o {output}"


rule get_nan_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=get_nan_labels,
    output:
        mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="nan",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule create_upsampled_coords_ref:
    input:
        seg=get_input_for_shape_inject,
    params:
        tight_crop_labels=lambda wildcards: config["tight_crop_labels"][wildcards.label],
        resample_res=lambda wildcards: config["laminar_coords_res"][wildcards.label],
        trim_padding="3mm",
    output:
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            space="corobl",
            label="{label}",
            hemi="{hemi}",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels {params.tight_crop_labels} -trim {params.trim_padding} -resample-mm {params.resample_res} -o {output}"


rule prep_dseg_for_laynii:
    input:
        dseg_tissue=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: " ".join(
            [
                str(lbl)
                for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                    "gm"
                ]
            ]
        ),
        src_labels=lambda wildcards: " ".join(
            [
                str(lbl)
                for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                    "src"
                ]
            ]
        ),
        sink_labels=lambda wildcards: " ".join(
            [
                str(lbl)
                for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                    "sink"
                ]
            ]
        ),
    output:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir,IO}",
            desc="laynii",
            label="{label}",
            space="corobl",
            hemi="{hemi}",
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -as DSEG -retain-labels {params.gm_labels} -binarize -scale 3 -popas GM -push DSEG -retain-labels {params.src_labels} -binarize -scale 2 -popas WM -push DSEG -retain-labels {params.sink_labels} -binarize -scale 1 -popas PIAL -push GM -push WM -add -push PIAL -add -o {output}"


rule laynii_layers_equidist:
    input:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir}",
            desc="laynii",
            label="{label}",
            space="corobl",
            hemi="{hemi}",
        ),
    output:
        equidist=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="equidist",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("laynii")
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN2_LAYERS  -rim dseg.nii.gz && "
        "cp dseg_metric_equidist.nii.gz {output.equidist}"


rule laynii_layers_equivol:
    input:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir}",
            desc="laynii",
            label="{label}",
            space="corobl",
            hemi="{hemi}",
        ),
    output:
        equivol=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="equivol",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("laynii")
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN2_LAYERS  -rim dseg.nii.gz -equivol && "
        "cp dseg_metric_equivol.nii.gz {output.equivol}"
