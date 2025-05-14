

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
        labelmap=get_input_for_shape_inject,
    params:
        labels=get_gm_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                desc="GM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    conda:
        conda_env("c3d")
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"



rule get_src_sink_mask:
    input:
        labelmap=get_input_for_shape_inject,
    params:
        labels=get_src_sink_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="{srcsink,src|sink}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
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
            root=root,
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
        sdt=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="sdt.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="{srcsink}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -sdt -o {output}"


rule get_nan_mask:
    input:
        labelmap=get_input_for_shape_inject,
    params:
        labels=get_nan_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="nan",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule prep_dseg_for_laynii:
    input:
        dseg_tissue=get_input_for_shape_inject,
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
        dseg_rim=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                dir="{dir,IO}",
                desc="laynii",
                label="{label}",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -as DSEG -retain-labels {params.gm_labels} -binarize -scale 3 -popas GM -push DSEG -retain-labels {params.src_labels} -binarize -scale 2 -popas WM -push DSEG -retain-labels {params.sink_labels} -binarize -scale 1 -popas PIAL -push GM -push WM -add -push PIAL -add -o {output}"


rule laynii_layers:
    input:
        dseg_rim=bids(
            root=root,
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
        equidist=temp(
            bids(
                root=root,
                datatype="coords",
                dir="{dir,IO}",
                label="{label}",
                suffix="coords.nii.gz",
                desc="equidist",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    shadow:
        "minimal"
    conda:
        conda_env("laynii")
    log:
        bids_log(
            "laynii_layers",
            **inputs.subj_wildcards,
            dir="{dir, IO}",
            label="{label}",
            hemi="{hemi}",
        ),
    group:
        "subj"
    shell:
        "cp {input} dseg.nii.gz && "
        "LN2_LAYERS  -rim dseg.nii.gz &> {log} && "
        # "LN2_LAYERS  -rim dseg.nii.gz -equivol &> {log} && "
        "cp dseg_metric_equidist.nii.gz {output.equidist}"
