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

rule morphclose_dg:
    input:
        dseg_tissue=get_labels_for_laplace,
    output:
        dseg_tissue=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="closeDG",
            space="corobl",
            hemi="{hemi}"
        ),
    group: 'subj'
    container:
        config["singularity"]["autotop"]
    shell: 
        "c3d {input} -as DSEG -retain-labels 8 -binarize -dilate 1 3x3x3vox -erode 1 3x3x3vox -scale 100 -push DSEG -max -replace 100 8 -o {output}"

rule prep_dseg_for_laynii_hipp:
    input:
        dseg_tissue=get_labels_for_laplace,
    params:
        gm_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in config["laplace_labels"][wildcards.dir]["gm_noDG"]]
        ),
        src_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in config["laplace_labels"][wildcards.dir]["src"]]
        ),
        sink_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in config["laplace_labels"][wildcards.dir]["sink"]]
        ),
    output:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir,IO}",
            desc="laynii",
            label="{autotop,hipp}",
            space="corobl",
            hemi="{hemi}"
        ),
    container:
        config["singularity"]["autotop"]
    group:
        "subj"
    shell:
        "c3d -background -1 {input} -as DSEG -retain-labels {params.gm_labels} -binarize -scale 3 -popas GM -push DSEG -retain-labels {params.src_labels} -binarize -scale 2 -popas WM -push DSEG -retain-labels {params.sink_labels} -binarize -scale 1 -popas PIAL -push GM -push WM -add -push PIAL -add -o {output}"


rule prep_dseg_for_laynii_dentate:
    input:
        dseg_tissue=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="closeDG",
            space="corobl",
            hemi="{hemi}"
        ),
    params:
        gm_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in [8]]
        ),
        src_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in [2,4,7,0] ]
        ),
        sink_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in [1]]
        ),
    output:
        dseg_rim=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            dir="{dir,IO}",
            desc="laynii",
            label="{autotop,dentate}",
            space="corobl",
            hemi="{hemi}"
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
            hemi="{hemi}"
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
            **inputs.subj_wildcards
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
            **inputs.subj_wildcards
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
    """Renzo implementation of equidist, using LN_GROW_LAYERS"""
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
            hemi="{hemi}"
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
            **inputs.subj_wildcards
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
        "cp dseg_metric_equidist.nii.gz {output.equidist}" # TODO: naming



rule laynii_equivol_renzo:
    """Renzo implementation of equivol, using LN_GROW_LAYERS"""
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
            hemi="{hemi}"
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
            **inputs.subj_wildcards
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
        "cp dseg_metric_equivol.nii.gz {output.equivol}" #-- TODO naming
