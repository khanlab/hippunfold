import itertools
import re

import pandas as pd

# first surface is used as a reference when generating average shape
ref_hemi = inputs[config["modality"]].zip_lists.get("hemi", ["L"])[0]


def get_cmd_templategen_subj_csv(wildcards, input, output):

    cmds = []

    # write to csv
    csv_line = ",".join(input.metric_nii)
    cmds.append(f"echo {csv_line} > {output.metrics_csv}")

    return " && ".join(cmds)


rule align_lr_unfold_2d:
    """we average L and R together, this ensures they are in the same space"""
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix="{metric}.nii.gz",
            space="unfold",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        img=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfoldLR",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c3d {input} {params.flip_per_hemi} -o {output}"


rule templategen_subj_csv:
    input:
        metric_nii=expand(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=config["new_atlas_metrics"],
            allow_missing=True,
        ),
    output:
        metrics_csv=temp(
            bids(
                root=root,
                datatype="anat",
                suffix="metrics.csv",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        cmd=get_cmd_templategen_subj_csv,
    shell:
        "{params.cmd}"


rule template_gen_combined_csv:
    input:
        metrics_csvs=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix="metrics.csv",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            **expand_hemi_atlas_gen(wildcards),
        ),
    output:
        metrics_csv=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
            **hemi_wildcard_atlas_gen(),
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        cmd=lambda wildcards, input, output: f"cat {input.metrics_csvs} > {output.metrics_csv}",
    shell:
        "{params.cmd}"


# antsMultivariateTemplateConstruction2.sh options and the set of metrics (mets),
# also the gen_atlas_reg_ants_tune wildcards
TUNE_KEYS = [
    "mets",
    "sim",
    "tx",
    "step",
    "updvar",
    "totvar",
    "tplstep",
    "iters",
    "q",
    "f",
    "s",
    "wt",
    "sharpen",
]
TUNE_WILDCARD_CONSTRAINTS = {key: "[a-zA-Z0-9.]+" for key in TUNE_KEYS}
TUNE_WILDCARD_CONSTRAINTS["mets"] = "[a-z]+"

# short codes for the metrics, concatenated for the mets wildcard (e.g. curvgyrthick)
METRIC_CODES = {
    "curvature": "curv",
    "gyrification": "gyr",
    "legacygyrification": "lgyr",
    "thickness": "thick",
    "myelin": "myelin",
}


def metrics_to_code(metrics):
    """mets wildcard value for a list of metrics, in the order of METRIC_CODES"""
    return "".join(code for metric, code in METRIC_CODES.items() if metric in metrics)


def code_to_metrics(code):
    """list of metrics from a mets wildcard value"""
    codes = re.findall(
        "|".join(sorted(METRIC_CODES.values(), key=len, reverse=True)), code
    )
    if "".join(codes) != code or len(set(codes)) != len(codes):
        raise ValueError(
            f"cannot parse metrics from {code}, use a combination of {list(METRIC_CODES.values())}"
        )
    return [m for m, c in METRIC_CODES.items() if c in codes]


def get_ants_template_params(overrides={}):
    """new_atlas_ants from the config and the new_atlas_metrics (as mets), with the
    default modality weights filled in"""
    params = {
        "mets": metrics_to_code(config["new_atlas_metrics"]),
        **config["new_atlas_ants"],
        **overrides,
    }
    if not params.get("wt"):
        params["wt"] = "x".join(["1"] * len(code_to_metrics(params["mets"])))
    return {key: str(params[key]) for key in TUNE_KEYS}


def get_ants_template_opts(params, num_modalities):
    """command-line options for antsMultivariateTemplateConstruction2.sh"""
    if len(params["wt"].split("x")) != num_modalities:
        raise ValueError(
            f"modality weights {params['wt']} must have {num_modalities} components"
        )
    if not (
        len(params["q"].split("x"))
        == len(params["f"].split("x"))
        == len(params["s"].split("x"))
    ):
        raise ValueError(
            f"q, f, s must have the same number of levels: {params['q']}, {params['f']}, {params['s']}"
        )
    sim = re.sub(r"^CC(\d+)$", r"CC[\1]", params["sim"])
    regularization = f"{params['step']},{params['updvar']},{params['totvar']}"
    if params["tx"] == "SyN":
        transform = f"SyN[{regularization}]"
    elif params["tx"] == "BSplineSyN":
        transform = f"BSplineSyN[{regularization},3]"
    else:
        raise ValueError(f"unsupported transform {params['tx']}, use SyN or BSplineSyN")
    return (
        f"-k {num_modalities} -m '{sim}' -t '{transform}' -g {params['tplstep']} "
        f"-i {params['iters']} -q {params['q']} -f {params['f']} -s {params['s']} "
        f"-w {params['wt']} -A {params['sharpen']}"
    )


def get_tune_path(suffix):
    return bids_atlas(
        root=root,
        template=config["new_atlas_name"],
        datatype="tune",
        label="{label}",
        **hemi_wildcard_atlas_gen(),
        **{key: f"{{{key}}}" for key in TUNE_KEYS},
        suffix=suffix,
    )


def get_tune_grid():
    """list of parameter dicts for the grid in new_atlas_tune: the product of the grid
    lists, with q/f/s varied together from levels, other keys from new_atlas_ants"""
    tune = config["new_atlas_tune"]
    grid = tune.get("grid", {}) or {}
    unknown = set(grid) - set(TUNE_KEYS)
    if unknown:
        raise ValueError(f"unknown keys in new_atlas_tune grid: {unknown}")
    if set(grid) & {"q", "f", "s"}:
        raise ValueError("set q, f and s in new_atlas_tune levels, not the grid")
    base = get_ants_template_params()
    # modality weights default to all 1 for each metric set, unless in the grid
    reset_wt = {} if "wt" in grid else {"wt": None}
    levels = tune.get("levels") or [{key: base[key] for key in ("q", "f", "s")}]
    keys = list(grid)
    combos = []
    for values in itertools.product(*[grid[key] for key in keys]):
        for level in levels:
            combos.append(
                get_ants_template_params(
                    {**reset_wt, **dict(zip(keys, values)), **level}
                )
            )
    return combos


def get_tune_files(suffix):
    hemis = (
        [{}]
        if config["new_atlas_subfields_from"] == "native"
        else [{"hemi": hemi} for hemi in config["hemi"]]
    )
    labels = config["new_atlas_tune"].get("labels") or config["autotop_labels"]
    return [
        get_tune_path(suffix).format(label=label, **hemi, **params)
        for label in labels
        for hemi in hemis
        for params in get_tune_grid()
    ]


def get_tune_report():
    return bids_atlas(
        root=root,
        template=config["new_atlas_name"],
        datatype="tune",
        desc="tune",
        suffix="report.html",
    )


rule gen_atlas_reg_ants:
    input:
        metrics_niis=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=config["new_atlas_metrics"],
            label=wildcards.label,
            **expand_hemi_atlas_gen(wildcards),
        ),
        metrics_csv=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            desc="combined",
            suffix="metrics.csv",
            **hemi_wildcard_atlas_gen(),
        ),
    output:
        avgtemplate_dir=directory(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                suffix="antstemplate",
                **hemi_wildcard_atlas_gen(),
            )
        ),
    threads: 4
    resources:
        mem_mb=scale_by_attempt(16000),
        runtime=scale_by_attempt(720),
    params:
        ants_opts=get_ants_template_opts(
            get_ants_template_params(), len(config["new_atlas_metrics"])
        ),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsMultivariateTemplateConstruction2.sh "
        " -d 2 -o {params.warp_prefix} -n 0 -l 0 -b 1 "
        " {params.ants_opts} {input.metrics_csv}"


def get_tune_metric_rows(wildcards):
    """one row of metric images per input (subject/hemi), for the metrics in the
    mets wildcard, in the same input order as template_gen_combined_csv"""
    per_metric = [
        inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix=f"{metric}.nii.gz",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            **expand_hemi_atlas_gen(wildcards),
        )
        for metric in code_to_metrics(wildcards.mets)
    ]
    return list(zip(*per_metric))


def get_tune_csv_path():
    return bids_atlas(
        root=root,
        template=config["new_atlas_name"],
        datatype="tune",
        label="{label}",
        **hemi_wildcard_atlas_gen(),
        mets="{mets}",
        desc="combined",
        suffix="metrics.csv",
    )


rule template_gen_combined_csv_tune:
    """same as template_gen_combined_csv, for the metrics in the mets wildcard"""
    input:
        metrics_niis=lambda wildcards: [
            f for row in get_tune_metric_rows(wildcards) for f in row
        ],
    output:
        metrics_csv=get_tune_csv_path(),
    wildcard_constraints:
        mets=TUNE_WILDCARD_CONSTRAINTS["mets"],
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        rows=get_tune_metric_rows,
    run:
        with open(output.metrics_csv, "w") as f:
            f.writelines(",".join(row) + "\n" for row in params.rows)


rule gen_atlas_reg_ants_tune:
    """same as gen_atlas_reg_ants, but with the antsMultivariateTemplateConstruction2.sh
    options and the set of metrics as wildcards (see new_atlas_ants in the config),
    for hyperparameter tuning
    """
    input:
        metrics_niis=lambda wildcards: [
            f for row in get_tune_metric_rows(wildcards) for f in row
        ],
        metrics_csv=get_tune_csv_path(),
    output:
        avgtemplate_dir=directory(get_tune_path("antstemplate")),
    wildcard_constraints:
        **TUNE_WILDCARD_CONSTRAINTS,
    threads: 4
    resources:
        mem_mb=scale_by_attempt(16000),
        runtime=scale_by_attempt(720),
    params:
        ants_opts=lambda wildcards: get_ants_template_opts(
            {key: getattr(wildcards, key) for key in TUNE_KEYS},
            len(code_to_metrics(wildcards.mets)),
        ),
        warp_prefix=lambda wildcards, output: f"{output.avgtemplate_dir}/",
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsMultivariateTemplateConstruction2.sh "
        " -d 2 -o {params.warp_prefix} -n 0 -l 0 -b 0 "
        " {params.ants_opts} {input.metrics_csv}"


def get_tune_subfields(wildcards):
    """subfield labels in the same order as the rows of the combined metrics csv
    (only the hipp has subfields)"""
    if wildcards.label != "hipp":
        return []
    return inputs[config["modality"]].expand(
        bids(
            root=root,
            datatype="anat",
            suffix="subfieldsfrom{subfields_from}.nii.gz".format(
                subfields_from=config["new_atlas_subfields_from"]
            ),
            space="unfoldLR2d",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        label=wildcards.label,
        **expand_hemi_atlas_gen(wildcards),
    )


rule eval_atlas_reg_ants_tune:
    """leave-one-out subfield dice, metric correlation and warp regularity
    for one template construction run"""
    input:
        avgtemplate_dir=get_tune_path("antstemplate"),
        metrics_csv=get_tune_csv_path(),
        subfields=get_tune_subfields,
    output:
        tsv=get_tune_path("eval.tsv"),
        png=get_tune_path("eval.png"),
    wildcard_constraints:
        **TUNE_WILDCARD_CONSTRAINTS,
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(30),
    params:
        tune_params=lambda wildcards: {
            "label": wildcards.label,
            **({"hemi": wildcards.hemi} if "hemi" in wildcards.keys() else {}),
            **{key: getattr(wildcards, key) for key in TUNE_KEYS},
        },
    script:
        "../scripts/eval_template_tuning.py"


rule tune_atlas_report:
    input:
        tsvs=lambda wildcards: get_tune_files("eval.tsv"),
        pngs=lambda wildcards: get_tune_files("eval.png"),
    output:
        tsv=get_tune_report().replace(".html", ".tsv"),
        html=get_tune_report(),
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(15),
    params:
        tune_keys=TUNE_KEYS,
        base_params=lambda wildcards: get_ants_template_params(),
    script:
        "../scripts/make_tuning_report.py"


rule copy_avgtemplate_warps:
    input:
        metric_nii=expand(
            bids(
                root=root,
                datatype="anat",
                suffix="{metric}.nii.gz",
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            metric=config["new_atlas_metrics"],
            allow_missing=True,
        ),
        avgtemplate_dir=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            suffix="antstemplate",
            **hemi_wildcard_atlas_gen(),
        ),
    output:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                suffix="invwarp.nii.gz",
                from_="unfold",
                to=config["new_atlas_name"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        glob_input_warp=lambda wildcards, input: "{avgtemplate_dir}/input*-{filename}-1Warp.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            filename=Path(input.metric_nii[0]).name.removesuffix(".nii.gz"),
        ),
        glob_input_invwarp=lambda wildcards, input: "{avgtemplate_dir}/input*-{filename}-1InverseWarp.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            filename=Path(input.metric_nii[0]).name.removesuffix(".nii.gz"),
        ),
    shell:
        "cp {params.glob_input_warp} {output.warp} && "
        "cp {params.glob_input_invwarp} {output.invwarp}"


rule unflip_avgtemplate_metric:
    """copy avgtemplate metric out of avgtemplate folder"""
    input:
        avgtemplate_dir=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            suffix="antstemplate",
            **hemi_wildcard_atlas_gen(),
        ),
    output:
        metric=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                desc="badhdr",
                suffix="{metric}.nii.gz",
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        in_metric=lambda wildcards, input: "{avgtemplate_dir}/template{i}.nii.gz".format(
            avgtemplate_dir=input.avgtemplate_dir,
            i=config["new_atlas_metrics"].index(wildcards.metric),
        ),
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c3d {params.in_metric} {params.flip_per_hemi} -o {output.metric}"


rule reset_header_2d_metric_nii:
    """adjusts header to match the original data
    (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="metricref.nii.gz",
        ),
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            desc="badhdr",
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric,[a-zA-Z0-9]+}.nii.gz",
        ),
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    script:
        "../scripts/set_metric_nii_header.py"


rule reset_header_2d_warp_atlasgen:
    """adjusts header to match the original data
    (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="metricref.nii.gz",
        ),
        nii=bids(
            root=root,
            datatype="warps",
            suffix="{warpsuffix}.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="warps",
            suffix="{warpsuffix,warp|invwarp}.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            desc="3D",
            **inputs.subj_wildcards,
        ),
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    script:
        "../scripts/set_metric_nii_header.py"


rule create_unfold_ref_2d:
    output:
        metric_ref=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                suffix="metricref.nii.gz",
            )
        ),
    shadow:
        "minimal"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        dims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["dims"][:2]
        )
        + "x1",
        voxdims=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["voxdims"]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["origin"][:2]
        )
        + "x"
        + str(
            float(config["unfold_vol_ref"][wildcards.label]["extent"][-1])
            * config["surf_thresholds"]["midthickness"]
            * -1.0
        ),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.label]["orient"],
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} {params.flip_per_hemi} -slice z 50%  -o {output}"


rule create_unfold_ref_2d_resampled:
    output:
        metric_ref=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                resample="{resample}",
                suffix="metricref.nii.gz",
            )
        ),
    shadow:
        "minimal"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        dims=lambda wildcards: "x".join(
            [
                str(int(float(d) * float(wildcards.resample) / 100.0))
                for d in config["unfold_vol_ref"][wildcards.label]["dims"][:2]
            ]
        ),
        voxdims=lambda wildcards: "x".join(
            [
                str(float(d) / (float(wildcards.resample) / 100.0))
                for d in config["unfold_vol_ref"][wildcards.label]["voxdims"][:2]
            ]
        ),
        origin=lambda wildcards: "x".join(
            config["unfold_vol_ref"][wildcards.label]["origin"][:2]
        ),
        orient=lambda wildcards: config["unfold_vol_ref"][wildcards.label]["orient"][:2],
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c2d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} {params.flip_per_hemi} -scale 0 -shift 1 -binarize  -o {output}"


rule gen_unfold_atlas_mesh:
    input:
        metric_ref=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            resample="{resample}",
            suffix="metricref.nii.gz",
        ),
    output:
        surf_gii=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                resample="{resample}",
                space="unfold",
                hemi="{hemi,R}",
                suffix="{surfname,midthickness|inner|outer}.surf.gii",
            )
        ),
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(15),
    params:
        z_level=get_unfold_z_level,
    script:
        "../scripts/gen_unfold_atlas_mesh.py"


rule gen_unfold_atlas_mesh_flip:
    input:
        surf_gii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            resample="{resample}",
            space="unfold",
            hemi="R",
            suffix="{surfname}.surf.gii",
        ),
    output:
        surf_gii=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                resample="{resample}",
                space="unfold",
                hemi="{hemi,L}",
                suffix="{surfname,midthickness|inner|outer}.surf.gii",
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        z_level=get_unfold_z_level,
    shell:
        "wb_command -surface-flip-lr {input} {output}"


def get_unfold_mesh_resample(wildcards):

    ind_resample = config["density_choices"].index(wildcards.density)
    resample = config["resample_factors"][ind_resample]

    return bids_atlas(
        root=root,
        template=config["new_atlas_name"],
        label="{label}",
        resample="{resample}",
        space="unfold",
        hemi="{hemi}",
        suffix="{surfname}.surf.gii",
    ).format(resample=resample, **wildcards)


rule update_unfold_mesh_metadata:
    input:
        get_unfold_mesh_resample,
    output:
        surf_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            space="unfold",
            hemi="{hemi}",
            suffix="{surfname,midthickness|inner|outer}.surf.gii",
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
        secondary_type=lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type="FLAT",
    shell:
        "wb_command -surface-flip-normals {input} {output} && "
        "wb_command -set-structure {output} {params.structure_type} -surface-type {params.surface_type}"
        " -surface-secondary-type {params.secondary_type}"


rule avgtemplate_metric_vol_to_surf:
    input:
        metric_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="{metric}.nii.gz",
        ),
        midthickness=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            den="{density}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            hemi="{hemi}",
            suffix="{metric}.shape.gii",
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
    shell:
        "wb_command -volume-to-surface-mapping {input.metric_nii} {input.midthickness} {output.metric_gii} -trilinear && "
        "wb_command -set-structure {output.metric_gii} {params.structure_type}"


rule warp_subj_unfold_surf_to_avg:
    input:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfold",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="unfold",
            to=config["new_atlas_name"],
            hemi="{hemi}",
            label="{label}",
            desc="3D",
            **inputs.subj_wildcards,
        ),
    output:
        surf_gii=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldavg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        cmd=get_cmd_warp_surface_2d_warp,
    shell:
        "{params.cmd}"


rule resample_subj_native_surf_to_avg:
    input:
        subj_native=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="corobl",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        subj_unfold=bids(
            root=root,
            datatype="surf",
            suffix="{surfname}.surf.gii",
            space="unfoldavg",
            den="native",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
        atlas_unfold=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            resample="{resample}",
            hemi="{hemi}",
            suffix="{surfname}.surf.gii",
        ),
    output:
        native_resampled=bids(
            root=root,
            datatype="surf",
            label="{label}",
            space="corobl",
            resample="{resample}",
            hemi="{hemi}",
            suffix="{surfname}.surf.gii",
            **inputs.subj_wildcards,
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    shell:
        "wb_command -surface-resample {input.subj_native} {input.subj_unfold} {input.atlas_unfold} BARYCENTRIC {output.native_resampled} -bypass-sphere-check"


rule vote_subfield_labels:
    """multi-label smoothing fusion: each subject's subfields (either defined from
    native (manual segs), or from the base unfolded atlas) are split into per-label
    channels, gaussian smoothed, warped to the average space with linear interp,
    averaged across subjects, and the label with highest average is chosen"""
    input:
        imgs=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="anat",
                suffix="subfieldsfrom{subfields_from}.nii.gz".format(
                    subfields_from=config["new_atlas_subfields_from"]
                ),
                space="unfoldLR2d",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            **expand_hemi_atlas_gen(wildcards),
        ),
        warps=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="warps",
                suffix="warp.nii.gz",
                from_="unfold",
                to=config["new_atlas_name"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            label=wildcards.label,
            **expand_hemi_atlas_gen(wildcards),
        ),
    output:
        subfields_voted=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                desc="subfields",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
                label="{label}",
            )
        ),
    resources:
        mem_mb=scale_by_attempt(4000),
        runtime=scale_by_attempt(30),
    params:
        sigma_mm=config.get("new_atlas_label_smoothing_sigma_mm", 0.3125),
    script:
        "../scripts/multilabel_smooth_voting.py"


rule reset_header_2d_subfields_nii:
    """adjusts header to match the original data
    (since this seems to get garbled in z by ants)"""
    input:
        ref_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            suffix="metricref.nii.gz",
        ),
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsunflip",
            suffix="dseg.nii.gz",
        ),
    output:
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsunflipfixhdr",
            suffix="dseg.nii.gz",
        ),
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    script:
        "../scripts/set_metric_nii_header.py"


rule unflip_avg_subfields_nii:
    input:
        nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfields",
            suffix="dseg.nii.gz",
        ),
    output:
        nii=temp(
            bids_atlas(
                root=root,
                template=config["new_atlas_name"],
                label="{label}",
                hemi="{hemi}",
                desc="subfieldsunflip",
                suffix="dseg.nii.gz",
            )
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        flip_per_hemi=lambda wildcards: config["unfold_vol_ref"][wildcards.label][
            "flip_per_hemi"
        ][wildcards.hemi],
    shell:
        "c3d {input} {params.flip_per_hemi} -o {output}"


rule import_avg_subfields_as_label:
    """read in volumetric subfield dseg, used for group_create_atlas only"""
    input:
        vol_dseg=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldsunflipfixhdr",
            suffix="dseg.nii.gz",
        ),
        label_list=Path(workflow.basedir)
        / "../resources/label_lut/labellist_withdg.txt",
    output:
        label_dseg=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            label="{label}",
            hemi="{hemi}",
            desc="subfieldswithlbl",
            suffix="dseg.nii.gz",
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    shell:
        "wb_command -volume-label-import {input.vol_dseg} {input.label_list} {output.label_dseg}"


rule avgtemplate_subfield_voted_vol_to_surf:
    input:
        subfields_nii=bids_atlas(
            root=root,
            template=config["new_atlas_name"],
            desc="subfieldswithlbl",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
            label="{label}",
        ),
        midthickness=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            space="unfold",
            den="{density}",
            hemi="{hemi}",
            suffix="midthickness.surf.gii",
        ),
    output:
        metric_gii=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            den="{density}",
            hemi="{hemi}",
            suffix="dseg.label.gii",
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        structure_type=lambda wildcards: get_structure(wildcards.hemi, wildcards.label),
    shell:
        "wb_command -volume-label-to-surface-mapping {input.subfields_nii} {input.midthickness} {output.metric_gii} && "
        "wb_command -set-structure {output.metric_gii} {params.structure_type}"


# input function for the rule aggregate
def get_atlas_inputs(wildcards):

    files = []

    for label in config["autotop_labels"]:
        for hemi in config["hemi"]:

            files.extend(
                expand(
                    bids_atlas(
                        root=get_atlas_dir(),
                        template=config["new_atlas_name"],
                        label=label,
                        hemi=hemi,
                        den="{density}",
                        suffix="{metric}.shape.gii",
                    ),
                    metric=config["new_atlas_metrics"],
                    density=config["density_choices"],
                )
            )
            files.extend(
                expand(
                    bids_atlas(
                        root=get_atlas_dir(),
                        template=config["new_atlas_name"],
                        label=label,
                        hemi=hemi,
                        den="{density}",
                        space="{space}",
                        suffix="{surfname}.surf.gii",
                    ),
                    space=["unfold", "native"],
                    surfname=["inner", "outer", "midthickness"],
                    density=config["density_choices"],
                )
            )
            if label == "hipp":
                files.extend(
                    expand(
                        bids_atlas(
                            root=get_atlas_dir(),
                            template=config["new_atlas_name"],
                            label="hipp",
                            hemi=hemi,
                            den="{density}",
                            suffix="dseg.label.gii",
                        ),
                        density=config["density_choices"],
                    )
                )

    return files


def density_to_resample(density_wildcard):
    ind_resample = config["density_choices"].index(density_wildcard)
    resample = config["resample_factors"][ind_resample]

    return resample


rule register_surf_to_ref:
    input:
        ref_surf=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                label="{label}",
                space="corobl",
                resample="{resample}".format(
                    resample=density_to_resample(wildcards.density)
                ),
                hemi="{hemi}",
                suffix="{surfname}.surf.gii",
                **inputs.subj_wildcards,
            ),
            surfname=wildcards.surfname,
            label=wildcards.label,
            hemi=wildcards.hemi,
        )[0],
        surf=lambda wildcards: bids(
            root=root,
            datatype="surf",
            label="{label}",
            space="corobl",
            resample="{resample}".format(
                resample=density_to_resample(wildcards.density)
            ),
            hemi=wildcards.hemi,
            suffix="{surfname}.surf.gii",
            **inputs.subj_wildcards,
        ),
    output:
        surf=bids(
            root=root,
            datatype="surf",
            label="{label}",
            space="coroblaligned",
            density="{density}",
            hemi="{hemi}",
            suffix="{surfname}.surf.gii",
            **inputs.subj_wildcards,
        ),
    shadow:
        "minimal"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    shell:
        "wb_command -surface-affine-regression {input.surf} {input.ref_surf} affine.xfm && "
        "wb_command -surface-apply-affine {input.surf} affine.xfm {output.surf}"


rule average_native_surfs:
    input:
        surfs=lambda wildcards: inputs[config["modality"]].expand(
            bids(
                root=root,
                datatype="surf",
                label="{label}",
                space="coroblaligned",
                density="{density}",
                hemi="{hemi}",
                suffix="{surfname}.surf.gii",
                **inputs.subj_wildcards,
            ),
            surfname=wildcards.surfname,
            label=wildcards.label,
            density=wildcards.density,
            hemi=ref_hemi,
        ),
    output:
        surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            hemi=ref_hemi,
            den="{density}",
            space="native",
            suffix="{surfname}.surf.gii",
        ),
    resources:
        mem_mb=scale_by_attempt(2000),
        runtime=scale_by_attempt(5),
    params:
        surf_args=lambda wildcards, input: " ".join(
            [f"-surf {surf}" for surf in input.surfs]
        ),
    shell:
        "wb_command -surface-average {output} {params.surf_args}"


rule flip_average_native_surf:
    input:
        surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            hemi=ref_hemi,
            den="{density}",
            space="native",
            suffix="{surfname}.surf.gii",
        ),
    output:
        surf=bids_atlas(
            root=get_atlas_dir(),
            template=config["new_atlas_name"],
            label="{label}",
            hemi="L" if ref_hemi == "R" else "R",
            den="{density}",
            space="native",
            suffix="{surfname}.surf.gii",
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        structure_type=lambda wildcards: get_structure(
            "L" if ref_hemi == "R" else "R", wildcards.label
        ),
    shell:
        "wb_command -surface-flip-lr {input} {output} && "
        "wb_command -set-structure {output} {params.structure_type}"


rule write_template_json:
    input:
        get_atlas_inputs,
    output:
        json=str(
            Path(
                bids(
                    root=get_atlas_dir(),
                    tpl=config["new_atlas_name"],
                )
            )
            / "template_description.json"
        ),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(5),
    params:
        template_description={
            "Identifier": config["new_atlas_name"],
            "metric_wildcards": config["new_atlas_metrics"],
            # tells unfoldreg that gyrification here is the smoothing-based index
            "gyrification_definition": "sgi",
            "label_wildcards": config["autotop_labels"],
            "hemi_wildcards": config["hemi"],
            "density_wildcards": config["density_choices"],
            "Authors": [""],
            "Acknowledgements": "",
            "BIDSVersion": "",
            "HowToAcknowledge": "",
            "License": "MIT",
            "Name": "HippUnfold surface-based template and subfield atlas",
            "RRID": "",
            "ReferencesAndLinks": [""],
            "TemplateFlowVersion": "",
        },
    script:
        "../scripts/write_template_json.py"


# TODO add rule to create zip file of new atlas and push to zenodo (need to add CLI to enable this)
