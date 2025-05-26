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


