# --- interactive surface viewers ---
#
# Self-contained HTML pages showing each surface in anatomical and unfolded
# space, coloured by its shape metrics and subfields (see lib/surface_viewer.py):
#   sub-*/qc/sub-*_den-{density}_desc-surfaces_viewer.html   (one subject)
#   group_den-{density}_desc-surfaces_viewer.html            (all subjects, group level)


def get_surface_viewer_records(wildcards, subject_level=True):
    """the surfaces, metrics and subfields of each subject, hemi and label"""
    component = inputs[config["modality"]]
    if subject_level:
        component = component.filter(**{k: wildcards[k] for k in inputs.subj_wildcards})
    # hemispheres per subject (and session): from the inputs if hemi is an
    # input wildcard (as for dsegtissue), otherwise config["hemi"]
    hemis = {}
    for values in zip(*component.zip_lists.values()):
        row = dict(zip(component.zip_lists, values))
        key = tuple((k, row[k]) for k in inputs.subj_wildcards)
        hemis.setdefault(key, [])
        for hemi in [row["hemi"]] if "hemi" in row else config["hemi"]:
            if hemi in config["hemi"] and hemi not in hemis[key]:
                hemis[key].append(hemi)
    records = []
    for key, subject_hemis in hemis.items():
        entities = dict(key)
        name = "_".join(
            f"{'sub' if k == 'subject' else 'ses' if k == 'session' else k}-{v}"
            for k, v in entities.items()
        )
        for hemi in sorted(subject_hemis):
            for label in config["autotop_labels"]:
                common = dict(
                    root=root, den=wildcards.density, hemi=hemi, label=label, **entities
                )
                records.append(
                    {
                        "subject": name,
                        "hemi": hemi,
                        "label": label,
                        "native": bids(
                            datatype="surf",
                            space="corobl",
                            suffix="midthickness.surf.gii",
                            **common,
                        ),
                        "unfold": bids(
                            datatype="surf",
                            space="unfold",
                            suffix="midthickness.surf.gii",
                            **common,
                        ),
                        "metrics": {
                            metric.split(".")[0]: bids(
                                datatype="metric", suffix=f"{metric}.gii", **common
                            )
                            for metric in get_gifti_metric_types(label)
                        },
                        "subfields": (
                            bids(
                                datatype="metric",
                                atlas=config["atlas"],
                                suffix="subfields.label.gii",
                                **common,
                            )
                            if label == "hipp"
                            else None
                        ),
                    }
                )
    return records


def surface_viewer_files(records):
    files = []
    for r in records:
        files += [r["native"], r["unfold"], *r["metrics"].values()]
        if r["subfields"]:
            files.append(r["subfields"])
    return files


rule surface_viewer_subject:
    """interactive viewer of one subject's surfaces and metrics"""
    input:
        lambda wildcards: surface_viewer_files(get_surface_viewer_records(wildcards)),
    output:
        html=bids(
            root=root,
            datatype="qc",
            den="{density}",
            desc="surfaces",
            suffix="viewer.html",
            **inputs.subj_wildcards,
        ),
    log:
        bids_log("surface_viewer_subject", den="{density}", **inputs.subj_wildcards),
    group:
        "subj"
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    params:
        records=get_surface_viewer_records,
        title="Surfaces",
        entities=lambda wildcards: {"den": wildcards.density},
    script:
        "../scripts/surface_viewer.py"


rule surface_viewer_group:
    """interactive viewer of every subject's surfaces and metrics, with a group
    overview of the unfolded maps for spotting outliers"""
    input:
        lambda wildcards: surface_viewer_files(
            get_surface_viewer_records(wildcards, subject_level=False)
        ),
    output:
        html=bids(
            root=root,
            prefix="group",
            den="{density}",
            desc="surfaces",
            suffix="viewer.html",
        ),
    log:
        bids_log("surface_viewer_group", den="{density}"),
    resources:
        mem_mb=scale_by_attempt(8000),
        runtime=scale_by_attempt(30),
    params:
        records=lambda wildcards: get_surface_viewer_records(
            wildcards, subject_level=False
        ),
        title="Group surfaces",
        entities=lambda wildcards: {"den": wildcards.density},
    script:
        "../scripts/surface_viewer.py"
