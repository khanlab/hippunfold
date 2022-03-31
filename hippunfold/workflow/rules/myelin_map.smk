
rule divide_t1_by_t2:
    input:
        t1=bids(
            root=work,
            datatype="anat",
            hemi="{hemi}",
            space="corobl",
            desc="preproc",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"]
        ),
        t2=bids(
            root=work,
            datatype="anat",
            hemi="{hemi}",
            space="corobl",
            desc="preproc",
            suffix="T2w.nii.gz",
            **config["subj_wildcards"]
        ),
    output:
        t1overt2=bids(
            root=work,
            datatype="anat",
            hemi="{hemi}",
            space="corobl",
            desc="preproc",
            suffix="T1wDividedByT2w.nii.gz",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input.t1} {input.t2} -divide -replace inf 1000 -inf -1000 NaN 0 -o {output}"


rule create_ribbon:
    input:
        ap_coord=bids(
            root=work,
            datatype="coords",
            dir="AP",
            label="{autotop}",
            suffix="coords.nii.gz",
            space="{space}",
            desc="laplace",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        ribbon=bids(
            root=work,
            datatype="anat",
            label="{autotop}",
            suffix="mask.nii.gz",
            space="{space}",
            desc="ribbon",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "c3d {input} -binarize -o {output}"


# sample on hipp & dg midthickness surfaces
rule sample_myelin_map_surf:
    """ samples myelin map on surf using corobl space """
    input:
        vol=bids(
            root=work,
            datatype="anat",
            space="corobl",
            hemi="{hemi}",
            desc="preproc",
            suffix="T1wDividedByT2w.nii.gz",
            **config["subj_wildcards"]
        ),
        mid=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="midthickness.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        inner=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="inner.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        outer=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="outer.surf.gii",
            space="corobl",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
        ribbon=bids(
            root=work,
            datatype="anat",
            label="{autotop}",
            suffix="mask.nii.gz",
            space="corobl",
            desc="ribbon",
            hemi="{hemi}",
            **config["subj_wildcards"]
        ),
    output:
        metric=bids(
            root=work,
            datatype="surf",
            den="{density}",
            suffix="myelin.shape.gii",
            space="{space}",
            hemi="{hemi}",
            label="{autotop}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -volume-to-surface-mapping {input.vol} {input.mid} {output.metric} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"