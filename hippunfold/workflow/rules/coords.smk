import os


def get_labels_for_laplace_hipp(wildcards):
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


def get_labels_for_laplace(wildcards):
    if wildcards.label == "dentate":
        seg = bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="combined",
            space="corobl",
            hemi="{hemi}",
            label="dentate",
        )
    else:

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


rule create_upsampled_coords_ref:
    input:
        seg=get_input_for_shape_inject,
    params:
        tight_crop_labels=lambda wildcards: config["tight_crop_labels"][wildcards.label],
        resample_res=lambda wildcards: config["laminar_coords_res"][wildcards.label],
    output:
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            space="corobl",
            hemi="{hemi}",
            label="{label}"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels {params.tight_crop_labels} -trim 10vox -resample-mm {params.resample_res} -o {output}"


rule dentate_template_dseg_to_PDsrcsink:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentatetissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentatePDsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels 1 2 -trim 2vox -split -foreach -resample-mm 0.05x0.05x0.05mm -smooth 2vox -endfor -merge  -o {output}"


rule dentate_template_dseg_to_APsrcsink:
    """ the AP labels for hipp are not adjacent to dentate boundaries, 
        so we extend those labels later -- so we label these APx instead"""
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentateAPxsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels 5 6 -trim 2vox -popas LBL "
        " -push LBL -retain-labels 5 -binarize -sdt -clip 1 inf -as SDT "
        "    -scale 0 -shift 1 -push SDT -divide -popas SRC_SCORE"
        " -push LBL -retain-labels 6 -binarize -sdt -clip 1 inf -as SDT "
        "    -scale 0 -shift 1 -push SDT -divide -popas SINK_SCORE "
        " -push LBL -scale 0 -push SRC_SCORE -push SINK_SCORE -vote "
        " -type uchar -replace 1 5 2 6 "
        " -int 0 -resample-mm 0.05x0.05x0.05mm -o {output}"
        #1 / SDT"
        #1 / SDT"


rule dentate_template_dseg_to_IOsrcsink:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentateIOsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels 1 -trim 2vox -split -foreach -resample-mm 0.05x0.05x0.05mm -smooth 2vox -endfor -merge -replace 1 3 0 4 -o {output}"


rule dentate_template_dseg_to_mask:
    input:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentatetissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        template_seg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentateIOmask",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -retain-labels 8 -binarize -trim 2vox -smooth 1vox  -resample-mm 0.05x0.05x0.05mm -smooth 1vox -threshold 0.25 inf 1 0   -o {output}"


rule inject_dentate_dseg:
    input:
        template_dseg=bids(
            root=work,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="dentate{dir}{dsegtype}",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}"
        ),
        matrix=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.txt",
            datatype="warps",
            desc="moments",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}",
        ),
        warp=bids(
            root=work,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        interp_opt="-ri NN",
    output:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="{autotop}{dir,APx|PD|IO}{dsegtype}",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.upsampled_ref} -rm {input.template_dseg} {output.dseg}  -r {input.warp} {input.matrix} "


rule expand_dentate_APsrcsink:
    """ the AP labels for hipp are not adjacent to dentate boundaries, 
        so we this expands them to the whole image """
    input:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentateAPxsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentateAPsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input}  -popas LBL "
        " -push LBL -retain-labels 5 -binarize -sdt -clip 1 inf -as SDT "
        "    -scale 0 -shift 1 -push SDT -divide -popas SRC_SCORE"
        " -push LBL -retain-labels 6 -binarize -sdt -clip 1 inf -as SDT "
        "    -scale 0 -shift 1 -push SDT -divide -popas SINK_SCORE "
        " -push LBL -scale 0 -push SRC_SCORE -push SINK_SCORE -vote "
        " -type uchar -replace 1 5 2 6 -o {output}"
        #1 / SDT"
        #1 / SDT"


rule upsample_dseg_to_dentate:
    input:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}"
        ),
        upsampled_ref=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="ref.nii.gz",
            desc="resampled",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}"
        ),
    params:
        interp_opt="-ri NN",
    output:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postprocresampled",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    threads: 8
    shell:
        "c3d -int 0 {input.upsampled_ref} {input.dseg} -reslice-identity -o {output.dseg}  "


rule combine_dentate_dseg:
    """ we have gm mask, and src/sink labels"""
    input:
        gm_mask=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentateIOmask",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
        srcsinkPD=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentatePDsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
        srcsinkIO=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentateIOsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
        srcsinkAP=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="dentateAPsrcsink",
            hemi="{hemi}",
            suffix="dseg.nii.gz"
        ),
    output:
        dseg=bids(
            root=work,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="combined",
            space="corobl",
            hemi="{hemi}",
            label="dentate",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/combine_dentate_dseg.py"


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
        "LN2_LAYERS  -rim dseg.nii.gz -equivol -iter_smooth 100 && "
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
        equidist=bids(
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
        "cp dseg_metric_equidist.nii.gz {output.equidist}"
        # TODO: naming


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
        "cp dseg_metric_equivol.nii.gz {output.equivol}"
        #-- TODO naming
