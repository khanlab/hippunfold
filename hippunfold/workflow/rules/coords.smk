import os


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
            hemi="{hemi}"
        ).format(**wildcards)
    return seg


def get_gm_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label]["IO"]["gm"]]
    )
    return lbl_list


def get_sink_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label][wildcards.dir]["sink"]]
    )
    return lbl_list


def get_src_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label][wildcards.dir]["src"]]
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
        gm_labels=get_gm_labels,
    output:
        mask=bids(
                root=work,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                desc="GM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -retain-labels {params} -binarize {output}"


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
                **inputs.subj_wildcards
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
                **inputs.subj_wildcards
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
                **inputs.subj_wildcards
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule equivolume_coords:
    input:
        source=bids(
                root=work,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="src",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
        ),
        gm=bids(
                root=work,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                desc="GM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
        ),
        edges=bids(
                root=work,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="nan",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
        ),
    params:
        script=os.path.join(workflow.basedir, "scripts/equivolume_coords.py"),
    output:
        coords=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            label="{label}",
            suffix="coords.nii.gz",
            desc="equivol",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        srcgm=bids(
                root=work,
                datatype="coords",
                suffix="mask.nii.gz",
                dir="{dir}",
                space="corobl",
                desc="srcGM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards
        ),
    group:
        "subj"
    resources:
        time=30,
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            dir="{dir}",
            hemi="{hemi}",
            label="{label}",
            suffix="equivolume.txt"
        ),
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.source} {input.gm} -add -o {output.srcgm} && "
        "c3d {output.srcgm} {input.edges} -add -o {output.srcgm} && "
        "python {params.script} {resources.tmpdir} {input.source} {output.srcgm} {output.coords} &> {log}"
