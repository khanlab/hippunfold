def get_labels_for_laplace(wildcards):
    if config["skip_inject_template_labels"]:
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
        ).format(**wildcards)
    return seg


def get_gm_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label]["IO"]["gm"]]
    )
    return lbl_list


def get_sink_labels(wildcards):
    lbl_list = " ".join(
        [
            str(lbl)
            for lbl in config["laplace_labels"][wildcards.label][wildcards.dir]["sink"]
        ]
    )
    return lbl_list


def get_src_labels(wildcards):
    lbl_list = " ".join(
        [
            str(lbl)
            for lbl in config["laplace_labels"][wildcards.label][wildcards.dir]["src"]
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


rule get_sink_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=get_sink_labels,
    output:
        mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="sink",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_src_mask:
    input:
        labelmap=get_labels_for_laplace,
    params:
        labels=get_src_labels,
    output:
        mask=bids(
            root=work,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="src",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


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
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


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
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -as DSEG -retain-labels {params.gm_labels} -binarize -scale 3 -popas GM -push DSEG -retain-labels {params.src_labels} -binarize -scale 2 -popas WM -push DSEG -retain-labels {params.sink_labels} -binarize -scale 1 -popas PIAL -push GM -push WM -add -push PIAL -add -o {output}"


rule laynii_layers:
    input:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir}",
            desc="laynii",
            label="{autotop}",
            space="corobl",
            hemi="{hemi}",
        ),
    output:
        equivol=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="equivol",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
        equidist=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{autotop}",
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
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN2_LAYERS  -rim dseg.nii.gz -equivol && "
        "cp dseg_metric_equidist.nii.gz {output.equidist} && "
        "cp dseg_metric_equivol.nii.gz {output.equivol}"


rule laynii_equidist_renzo:
    """Renzo implementation of equidist, using LN_GROW_LAYERS.  TODO: fix file names"""
    input:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir}",
            desc="laynii",
            label="{autotop}",
            space="corobl",
            hemi="{hemi}",
        ),
    output:
        equivol=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="equidistrenzo",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN_GROW_LAYERS  -rim dseg.nii.gz && "
        "cp dseg_metric_equidist.nii.gz {output.equidist}"


rule laynii_equivol_renzo:
    """Renzo implementation of equivol, using LN_GROW_LAYERS. TODO: fix filenames"""
    input:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir}",
            desc="laynii",
            label="{autotop}",
            space="corobl",
            hemi="{hemi}",
        ),
    output:
        equivol=bids(
            root=work,
            datatype="coords",
            dir="{dir,IO}",
            label="{autotop}",
            suffix="coords.nii.gz",
            desc="equivolrenzo",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN_GROW_LAYERS  -rim dseg.nii.gz -N 1000 -vinc 60 -threeD && "
        "LN_LEAKY_LAYERS  -rim dseg.nii.gz -nr_layers 1000 -iterations 100 && "
        "LN_LOITUMA  -equidist sc_rim_layers.nii -leaky sc_rim_leaky_layers.nii -FWHM 1 -nr_layers 10 && "
        "cp dseg_metric_equivol.nii.gz {output.equivol}"
