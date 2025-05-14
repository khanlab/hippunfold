# populate the HIPPUNFOLD_CACHE_DIR folder as needed
from lib import utils as utils

download_dir = utils.get_download_dir()

def get_inputs(): # TODO: swap this in
    if config["modality"] == "multires":
        return inputs["T1w"] + inputs["T2w"] + inputs["FLAIR"]
    else:
        return inputs[config["modality"]]   


rule import_any_modality:
    input:
        inputs[config["modality"]].path,
    output:
        bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        )
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule lamareg_to_template:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["T1w"].format(**wildcards),
    output:
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invaffine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
    shadow:
        "minimal"
    group:
        "subj" 
    log:
        bids_log(
            "lamareg_to_template",
            **inputs.wildcards,
        ),
    shell:
        "lamar generate-warpfield --ants-threads 4 --synthseg-thread 4 --fixed {params.ref} --moving {input.img} --affine {output.affine} --inverse-affine {output.invaffine} --warpfield {output.warp} --inverse-warpfield {output.invwarp} --inverse-output-parc tmp0.nii.gz --moving-parc tmp1.nii.gz --fixed-parc tmp2.nii.gz --registered-parc tmp3.nii.gz --output-parc tmp4.nii.gz &> {log}" # these SHOULD be removed in lamareg soon


rule apply_transforms:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_dir=Path(download_dir) / "template" / config["template"],
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
    params:
        ref=lambda wildcards, input: Path(input.template_dir)
        / config["template_files"][config["template"]]["Mask_crop"].format(**wildcards),
    output:
        img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    group:
        "subj"
    log:
        bids_log(
            "apply_transforms",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    shell:
        "lamar apply-warp --affine {input.affine} --warp {input.warp} --moving {input.img} --reference {params.ref} --output {output.img} &> {log}"


# TODO: refine registrations using the raw images and the above initializations


rule template_xfm_itk2ras:
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            **inputs.subj_wildcards,
            suffix="xfm.mat",
            from_="{modality}",
            to="corobl",
            type_="itk",
        ),
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                **inputs.subj_wildcards,
                suffix="xfm.mat",
                from_="{modality}",
                to="corobl",
                type_="ras",
            )
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -o {output}"


rule superres_inputs:
    input:
        bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    output:
        bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix="preproc.nii.gz",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        """
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            cp {input} {output}
        else
            c3d {input} -add -o {output}
        fi
        """