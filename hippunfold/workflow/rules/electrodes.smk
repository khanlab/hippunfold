rule import_electrodes_tsv:
    input:
        in_tsv=partial(get_single_bids_input, component="electrodes"),
    output:
        bids(
                root=root,
                datatype="ieeg",
                **inputs.subj_wildcards,
                suffix="electrodes.tsv",
            )
    group:
        "subj"
    shell:
        "cp {input} {output}"

 
rule import_electrodes_ref:
    input:
        in_img=partial(get_single_bids_input, component="electrodes_ref"),
    output:
        temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="electroderef.nii.gz",
            )
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"

rule reg_modality_to_electroderef:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="{modality}.nii.gz",
            desc="preproc",
        ),
        ref=bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="electroderef.nii.gz")
        
    output:
        warped=temp(bids(
            root=root,
            datatype="anat",
            suffix="{modality}.nii.gz",
            desc="preproc",
            space="electroderef",
            **inputs.subj_wildcards,
        )),
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.txt",
                from_="{modality}",
                to="electroderef",
                desc="rigid",
                type_="ras",
            )
        ),
    log:
        bids_log(
            "reg_modality_to_electroderef",
            modality='{modality}',
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["autotop"]
    conda:
        conda_env("niftyreg")
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac &> {log}"

rule annotate_electrodes_table:
    input:
        tsv=bids(
            root=root,
            datatype="ieeg",
            suffix="electrodes.tsv",
            **inputs.subj_wildcards,
        ),
        xfm=bids(
            root=root,
            datatype="warps",
            suffix="xfm.txt",
            from_=config["modality"],
            to="electroderef",
            desc="rigid",
            type_="ras",
            **inputs.subj_wildcards,
        ),
        surfs=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="{surfname}.surf.gii",
                space="{space}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            # surfname=["midthickness", "inner", "outer"],
            surfname=["midthickness"],
            hemi=["L", "R"],
            label=config['autotop_labels'],
            space=[config["modality"]],
            allow_missing=True,
        ),
        label_left=bids(
            root=root,
            datatype="metric",
            den="{density}",
            atlas=config["atlas"],
            suffix="subfields.label.gii",
            hemi="L",
            label="hipp",
            **inputs.subj_wildcards,
        ),
        label_right=bids(
            root=root,
            datatype="metric",
            den="{density}",
            atlas=config["atlas"],
            suffix="subfields.label.gii",
            hemi="R",
            label="hipp",
            **inputs.subj_wildcards,
        ),
    output:
        annotated=bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            desc="annotated",
            suffix="electrodes.tsv",
            **inputs.subj_wildcards,
        ),
    log:
        bids_log("annotate_electrodes_table",
                den="{density}",
                **inputs.subj_wildcards),
    conda:
        conda_env("pyvista"),
    container:
        config["singularity"]["autotop"],
    group:
        "subj",
    script:
        "../scripts/annotate_electrodes_table.py"


rule index_electrode_vertex_hits:
    input:
        annotated=bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            desc="annotated",
            suffix="electrodes.tsv",
            **inputs.subj_wildcards,
        ),
    output:
        hits=bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            suffix="electrodehits.json",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("pyunfold"),
    script:
        "../scripts/index_electrode_vertex_hits.py"

rule write_surface_label_gii:
    input:
        hits=bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            suffix="electrodehits.json",
            **inputs.subj_wildcards,
        ),
        surfs=expand(
            bids(
                root=root,
                datatype="surf",
                den="{density}",
                suffix="midthickness.surf.gii",
                space=config["modality"],
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=["L", "R"],
            label=config['autotop_labels'],
            allow_missing=True
        ),
    output:
        labelgii=expand(
            bids(
                root=root,
                datatype="ieeg",
                den="{density}",
                suffix="electrodes.label.gii",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            ),
            hemi=["L", "R"],
            label=config['autotop_labels'],
            allow_missing=True
        ),
    conda:
        conda_env("pyvista"),
    script:
        "../scripts/write_electrodes_label_gii.py"


rule electrode_label_gii_to_roi:
    input:
        bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            suffix="electrodes.label.gii",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        bids(
            root=root,
            datatype="ieeg",
            den="{density}",
            suffix="electrodes.shape.gii",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    conda:
         conda_env("workbench")
    group:
        "subj"
    shell:
        'wb_command -gifti-all-labels-to-rois {input} 1 {output}'


#should aggregate over  subjects using  metric-merge

