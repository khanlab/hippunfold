import os


def get_cmd_laplace_coords():
    if config["skip_inject_template_labels"]:
        cmd = "../scripts/laplace_coords.py"
    else:
        cmd = "../scripts/laplace_coords_withinit.py"
    return cmd


def get_labels_for_laplace(wildcards):
    if config["skip_inject_template_labels"]:
        seg = get_input_for_shape_inject(wildcards)
    else:
        seg = bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
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
                **config["subj_wildcards"],
                dir="{dir}",
                label="hipp",
                suffix="coords.nii.gz",
                desc="init",
                space="corobl",
                hemi="{hemi}"
            ),
        )
    return files


rule laplace_coords_hipp:
    input:
        unpack(get_inputs_laplace),
    params:
        gm_labels=lambda wildcards: config["laplace_labels"][wildcards.dir]["gm"],
        src_labels=lambda wildcards: config["laplace_labels"][wildcards.dir]["src"],
        sink_labels=lambda wildcards: config["laplace_labels"][wildcards.dir]["sink"],
        convergence_threshold=1e-5,
        max_iters=10000,
    output:
        coords=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            label="hipp",
            suffix="coords.nii.gz",
            desc="laplace",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    resources:
        time=30,
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            dir="{dir}",
            hemi="{hemi}",
            suffix="laplace-hipp.txt"
        ),
    container:
        config["singularity"]["autotop"]
    script:
        get_cmd_laplace_coords()


rule laplace_coords_dentate:
    input:
        coords=bids(
            root=work,
            datatype="coords",
            **config["subj_wildcards"],
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
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    resources:
        time=30,
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
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
            **config["subj_wildcards"]
        ),
        innerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="SRLM",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
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
            **config["subj_wildcards"]
        ),
        innerbin=bids(
            root=work,
            datatype="coords",
            dir="{dir}",
            desc="SRLM",
            suffix="mask.nii.gz",
            space="corobl",
            hemi="{hemi}",
            **config["subj_wildcards"]
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
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    resources:
        time=30,
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            dir="{dir}",
            hemi="{hemi}",
            suffix="equivolume.txt"
        ),
    container:
        config["singularity"]["autotop"]
    shell:
        "python {params.script} {resources.tmpdir} {input.innerbin} {input.outerbin} {output.coords} &> {log}"
