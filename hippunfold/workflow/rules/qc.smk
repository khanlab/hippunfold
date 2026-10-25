rule qc_reg_to_template:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="{modality}.nii.gz",
            space=config["template"],
            desc="affine",
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                **inputs.subj_wildcards,
                suffix="regqc.png",
                from_="{modality}",
                to=config["template"],
            ),
            caption="../report/t1w_template_regqc.rst",
            category="Registration QC",
        ),
    group:
        "preproc"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    params:
        ref=lambda wildcards, input: str(
            Path(input.template_dir)
            / config["template_files"][config["template"]][wildcards.modality].format(
                **wildcards
            )
        ),
    script:
        "../scripts/vis_regqc.py"


rule get_subfield_vols_subj:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=lambda wildcards: inputs[config["modality"]]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                **inputs.subj_wildcards,
                datatype="anat",
                hemi="{hemi}",
                space="{crop_ref_spaces}",
                desc="subfields",
                atlas="{atlas}",
                label="hipp",
                suffix="dseg.nii.gz",
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
        lookup_tsv=Path(workflow.basedir) / "../resources/label_lut/subfields_dseg.tsv",
    output:
        tsv=bids(
            root=root,
            datatype="anat",
            space="{crop_ref_spaces}",
            desc="subfields",
            atlas="{atlas}",
            suffix="volumes.tsv",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(1),
    script:
        "../scripts/gen_volume_tsv.py"


rule plot_subj_subfields:
    input:
        tsv=bids(
            root=root,
            datatype="anat",
            space="{crop_ref_spaces}",
            desc="subfields",
            atlas="{atlas}",
            suffix="volumes.tsv",
            **inputs.subj_wildcards,
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                space="{crop_ref_spaces}",
                desc="subfields",
                atlas="{atlas}",
                suffix="volumes.png",
                **inputs.subj_wildcards,
            ),
            caption="../report/subj_volume_plot.rst",
            category="Subfield Volumes",
        ),
    group:
        "subj"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(1),
    script:
        "../scripts/plot_subj_subfields.py"


def get_bg_img_for_subfield_qc(wildcards):
    if config["modality"] == "hippb500":
        return bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="hippb500.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        )
    elif config["modality"] == "dsegtissue":
        # blank image as bg
        return bids(
            root=root,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        )

    elif config["modality"][:3] == "seg":
        bg_modality = config["modality"][3:]
        return bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{bg_modality}.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        )

    else:
        bg_modality = config["modality"]
        return bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{bg_modality}.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **inputs.subj_wildcards,
        )


rule qc_subfield:
    input:
        img=get_bg_img_for_subfield_qc,
        seg=bids(
            root=root,
            datatype="anat",
            suffix="dseg.nii.gz",
            desc="subfields",
            space="{space}",
            hemi="{hemi}",
            atlas="{atlas}",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="dseg.png",
                desc="subfields",
                space="{space}",
                hemi="{hemi}",
                atlas="{atlas}",
                **inputs.subj_wildcards,
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(2),
    script:
        "../scripts/vis_qc_dseg.py"


rule qc_subfield_surf:
    input:
        surf=bids(
            root=root,
            datatype="surf",
            suffix="midthickness.surf.gii",
            den="{density}",
            space="{ref_spaces}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="midthickness.surf.png",
                den="{density}",
                desc="subfields",
                space="{ref_spaces}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "surf"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(2),
    script:
        "../scripts/vis_qc_surf.py"


rule concat_subj_vols_tsv:
    """Concatenate all subject tsv files into a single tsv"""
    input:
        tsv=lambda wildcards: inputs[get_modality_suffix(config["modality"])].expand(
            bids(
                root=root,
                datatype="anat",
                desc="subfields",
                space="{space}",
                atlas="{atlas}",
                suffix="volumes.tsv",
                **inputs.subj_wildcards,
            ),
            subject=inputs.subjects,
            session=inputs.sessions,
            space=wildcards.space,
            atlas=wildcards.atlas,
        ),
    output:
        tsv=bids(
            root=root,
            prefix="group",
            desc="subfields",
            space="{space}",
            atlas="{atlas}",
            from_="{modality}",
            suffix="volumes.tsv",
        ),
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    script:
        "../scripts/concat_tsv.py"


# --- QC checks and reports ---
#
# Each QC step (e.g. srcsink) writes QC metrics JSON files (see lib/qc.py),
# usually from the rule that computes what is being checked. Those are then
# visualized by a report rule, producing:
#   sub-*/qc/sub-*_desc-{qc}_qc.html
# and verified by qc_verify, which also requires the report so that it is
# always available to inspect when a check fails:
#   sub-*/qc/sub-*_desc-{qc}_qcchecks.tsv
# To add a QC step: write its metrics JSONs, add a report rule, and register
# a function returning its metrics JSONs in qc_metrics_inputs.


def get_srcsink_qc_metrics(wildcards):
    return (
        inputs[config["modality"]]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                datatype="qc",
                suffix="qcmetrics.json",
                desc="srcsink",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            label=config["autotop_labels"],
        )
    )


def get_srcsink_qc_surf_inputs(wildcards, **kwargs):
    return (
        inputs[config["modality"]]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                hemi="{hemi}",
                label="{label}",
                den="native",
                **kwargs,
                **inputs.subj_wildcards,
            ),
            hemi=config["hemi"],
            label=config["autotop_labels"],
        )
    )


qc_metrics_inputs = {
    "srcsink": get_srcsink_qc_metrics,
}


rule qc_report_srcsink:
    """HTML report of the AP/PD src/sink terminals on the boundary of each surface"""
    input:
        metrics=get_srcsink_qc_metrics,
        surf=lambda wildcards: get_srcsink_qc_surf_inputs(
            wildcards,
            datatype="surf",
            space="corobl",
            suffix="midthickness.surf.gii",
        ),
        coords_ap=lambda wildcards: get_srcsink_qc_surf_inputs(
            wildcards,
            datatype="metric",
            dir="AP",
            desc="laplace",
            suffix="coords.shape.gii",
        ),
        coords_pd=lambda wildcards: get_srcsink_qc_surf_inputs(
            wildcards,
            datatype="metric",
            dir="PD",
            desc="laplace",
            suffix="coords.shape.gii",
        ),
    output:
        html=bids(
            root=root,
            datatype="qc",
            desc="srcsink",
            suffix="qc.html",
            **inputs.subj_wildcards,
        ),
    log:
        bids_log("qc_report_srcsink", **inputs.subj_wildcards),
    group:
        "subj"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    script:
        "../scripts/qc_report_srcsink.py"


rule qc_verify:
    """Fails if any error-severity QC check failed (after the report is made)"""
    input:
        metrics=lambda wildcards: qc_metrics_inputs[wildcards.qc](
            {k: v for k, v in wildcards.items() if k != "qc"}
        ),
        report=bids(
            root=root,
            datatype="qc",
            desc="{qc}",
            suffix="qc.html",
            **inputs.subj_wildcards,
        ),
    output:
        tsv=bids(
            root=root,
            datatype="qc",
            desc="{qc}",
            suffix="qcchecks.tsv",
            **inputs.subj_wildcards,
        ),
    log:
        bids_log("qc_verify", desc="{qc}", **inputs.subj_wildcards),
    wildcard_constraints:
        qc="|".join(qc_metrics_inputs),
    group:
        "subj"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        warn_only=config.get("qc_warn_only", False),
    script:
        "../scripts/qc_verify.py"
