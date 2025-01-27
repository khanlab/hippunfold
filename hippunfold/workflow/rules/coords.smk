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
                hemi="{hemi}"
            ),
        )
    return files


rule laplace_coords_dentate:
    input:
        coords=bids(
            root=work,
            datatype="coords",
            **inputs.subj_wildcards,
            dir="{dir}",
            label="dentate",
            suffix="coords.nii.gz",
            desc="init",
            space="corobl",
            hemi="{hemi}"
        ),
    output:
        coords=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            label="dentate",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
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
            suffix="laplace-dentate.txt"
        ),
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output}"


rule prep_equivolume_coords:
    input:
        get_labels_for_laplace,
    params:
        src_labels=lambda wildcards: config["laplace_labels"][wildcards.dir]["src"],
    output:
        outerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="all",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        innerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="SRLM",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    log:
        bids(
            root="logs",
            **inputs.subj_wildcards,
            dir="{dir}",
            hemi="{hemi}",
            suffix="binarize.txt"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/prep_equivolume_coords.py"


rule equivolume_coords:
    input:
        outerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="all",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
        innerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="SRLM",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **inputs.subj_wildcards
        ),
    params:
        script=os.path.join(workflow.basedir, "scripts/equivolume_coords.py"),
    output:
        coords=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            label="hipp",
            suffix="coords.nii.gz",
            desc="equivol",
            space="corobl",
            hemi="{hemi}",
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
            suffix="equivolume.txt"
        ),
    container:
        config["singularity"]["autotop"]
    shell:
        "python {params.script} {resources.tmpdir} {input.innerbin} {input.outerbin} {output.coords} &> {log}"
