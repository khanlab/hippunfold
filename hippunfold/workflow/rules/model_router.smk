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
    
    return config["resource_urls"]["nnunet_model_version"].get(model_name, "v1")


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
    conda:
        "../envs/curl.yaml"
    shell:
        "mkdir -p {params.model_dir} && curl -L https://{params.url} -o {output}"


def get_inference_output_for_routing(wildcards):
    """Return appropriate input based on model version from config"""
    version = get_model_version()
    
    bids_path = bids(
        root=root,
        datatype="anat",
        **inputs.subj_wildcards,
        suffix="dseg.nii.gz",
        desc="nnunetv2" if version == "v2" else "nnunetv1",
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
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="nnunet",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    shell:
        "cp {input.seg} {output.nnunet_seg}"
