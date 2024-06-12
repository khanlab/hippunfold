rule qc_reg_to_template:
    input:
        flo=bids(
            root=work,
            datatype="anat",
            **config["subj_wildcards"],
            suffix="{native_modality}.nii.gz",
            space=config["template"],
            desc="affine"
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: str(
            Path(input.template_dir)
            / config["template_files"][config["template"]][
                wildcards.native_modality
            ].format(**wildcards)
        ),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                **config["subj_wildcards"],
                suffix="regqc.png",
                from_="{native_modality}",
                to=config["template"]
            ),
            caption="../report/t1w_template_regqc.rst",
            category="Registration QC",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/vis_regqc.py"


rule get_subfield_vols_subj:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            bids(
                root=root,
                **config["subj_wildcards"],
                datatype="anat",
                hemi="{hemi}",
                space="{crop_ref_spaces}",
                desc="subfields",
                atlas="{atlas}",
                suffix="dseg.nii.gz"
            ),
            hemi=config["hemi"],
            allow_missing=True,
        ),
        atlas_dir=lambda wildcards: Path(download_dir) / "atlas" / wildcards.atlas,
    params:
        lookup_tsv=lambda wildcards, input: Path(input.atlas_dir)
        / config["atlas_files"][wildcards.atlas]["lut"].format(**wildcards),
    group:
        "subj"
    output:
        tsv=bids(
            root=root,
            datatype="anat",
            space="{crop_ref_spaces}",
            desc="subfields",
            atlas="{atlas}",
            suffix="volumes.tsv",
            **config["subj_wildcards"]
        ),
    container:
        config["singularity"]["autotop"]
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
            **config["subj_wildcards"]
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
                **config["subj_wildcards"]
            ),
            caption="../report/subj_volume_plot.rst",
            category="Subfield Volumes",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/plot_subj_subfields.py"


def get_bg_img_for_subfield_qc(wildcards):
    if config["modality"] == "hippb500":
        return bids(
            root=work,
            datatype="anat",
            desc="preproc",
            suffix="hippb500.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **config["subj_wildcards"],
        )
    elif config["modality"] == "cropseg":
        # blank image as bg
        return bids(
            root=work,
            datatype="warps",
            suffix="cropref.nii.gz",
            space="{space}",
            hemi="{hemi}",
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"],
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
            **config["subj_wildcards"]
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
                **config["subj_wildcards"]
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
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
            label="{autotop}",
            **config["subj_wildcards"]
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
                label="{autotop}",
                **config["subj_wildcards"]
            ),
            caption="../report/subfield_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/vis_qc_surf.py"


rule concat_subj_vols_tsv:
    """Concatenate all subject tsv files into a single tsv"""
    input:
        tsv=lambda wildcards: expand(
            bids(
                root=root,
                datatype="anat",
                desc="subfields",
                space="{space}",
                atlas="{atlas}",
                suffix="volumes.tsv",
                **config["subj_wildcards"]
            ),
            subject=config["input_lists"][get_modality_key(config["modality"])][
                "subject"
            ],
            session=config["sessions"],
            space=wildcards.space,
            atlas=wildcards.atlas,
        ),
    group:
        "aggregate"
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
    container:
        config["singularity"]["autotop"]
    script:
        "../scripts/concat_tsv.py"
