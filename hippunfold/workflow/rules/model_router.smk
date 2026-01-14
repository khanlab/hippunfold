"""
Model routing rules for nnUNet v1 and v2
Routes inference to the appropriate nnUNet version based on config
"""
from pathlib import Path


def get_model_version():
    """Get nnUNet version from config file"""
    if config["force_nnunet_model"]:
        model_name = config["force_nnunet_model"]
    else:
        model_name = config["modality"]
    



    version_map = config["resource_urls"].get("nnunet_model_version", {})
    return version_map.get(model_name, "v1")


def get_model_tar():
    """Get path to the model tar file"""
    if config["force_nnunet_model"]:
        model_name = config["force_nnunet_model"]
    else:
        model_name = config["modality"]

    local_tar = config["resource_urls"]["nnunet_model"].get(model_name, None)
    if local_tar == None:
        print(f"ERROR: {model_name} does not exist in nnunet_model in the config file")

    return str((Path(download_dir) / "model" / Path(local_tar).name).absolute())


rule download_nnunet_model:
    """Download nnUNet model from remote URL"""
    params:
        url=(
            config["resource_urls"]["nnunet_model"][config["force_nnunet_model"]]
            if config["force_nnunet_model"]
            else config["resource_urls"]["nnunet_model"][config["modality"]]
        ),
        model_dir=Path(download_dir) / "model",
    output:
        model_tar=get_model_tar(),
    shell:
        "mkdir -p {params.model_dir} && curl -L https://{params.url} -o {output}"


def get_inference_output_for_routing(wildcards):
    """Return appropriate input based on model version from config"""
    version = get_model_version()
    
    if version == "synth_v1":
        desc = "synthseg"
    elif version == "v2":
        desc = "nnunetv2"
    else:
        desc = "nnunetv1"  # Changed to nnunetv1 to match new nnunet.smk
    
    bids_path = bids(
        root=root,
        datatype="anat",
        **config["subj_wildcards"],
        suffix="dseg.nii.gz",
        desc=desc,
        space="corobl",
        hemi="{hemi}",
    )
    
    return bids_path.format(**wildcards)


rule route_nnunet_inference:
    """Route to the appropriate nnUNet version based on config"""
    input:
        seg=get_inference_output_for_routing,
        model_tar=get_model_tar(),  # Ensure model is downloaded
    output:
        nnunet_seg=temp(
            bids(
                root=work,
                datatype="anat",
                **config["subj_wildcards"], # Changed from inputs.subj_wildcards
                suffix="dseg.nii.gz",
                desc="tissues",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    shell:
        "cp {input.seg} {output.nnunet_seg}"
