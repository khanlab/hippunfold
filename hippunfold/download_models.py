#!/usr/bin/env python3
import argparse
import errno
import os

import requests
import yaml
from appdirs import AppDirs


def get_model_dict():
    # get list of model urls
    model_cfg = os.path.join(
        os.path.dirname(__file__), "config", "nnunet_model_urls.yml"
    )
    with open(model_cfg, "r") as cfg:
        model_dict = yaml.load(cfg, Loader=yaml.FullLoader)
    return model_dict


def parse_args(model_dict):
    parser = argparse.ArgumentParser(
        prog="hippunfold_download_models",
        description="Tool for downloading U-net models for hippunfold",
    )

    parser.add_argument("--models", nargs="+", dest="models", choices=model_dict.keys())
    args = parser.parse_args()
    return args


def main():

    # get the model dict first, so we know what to parse
    model_dict = get_model_dict()
    inputs = parse_args(model_dict)

    if "HIPPUNFOLD_CACHE_DIR" in os.environ.keys():
        print(
            f"HIPPUNFOLD_CACHE_DIR defined, using: {os.environ['HIPPUNFOLD_CACHE_DIR']}"
        )
        download_dir = os.environ["HIPPUNFOLD_CACHE_DIR"]
    else:
        print(f"HIPPUNFOLD_CACHE_DIR not defined, using default location")
        # create local download dir if it doesn't exist
        dirs = AppDirs("hippunfold", "khanlab")
        download_dir = dirs.user_cache_dir

    try:
        os.mkdir(download_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    if inputs.models == None:
        models = model_dict.keys()
    else:
        models = inputs.models

    for modality in models:
        url = model_dict[modality]
        tarfile = url.split("/")[-1]
        local_path = os.path.join(download_dir, tarfile)

        # add ?dl=1 to url
        url = "".join([url, "?dl=1"])

        # if it doesn't exist, download the file
        if not os.path.exists(local_path):
            # download it:
            print(f"Downloading {modality} model...")
            print(f"   url = {url}")
            print(f"   dest = {local_path}")
            r = requests.get(url, allow_redirects=True, stream=True)
            with open(local_path, "wb") as f:
                f.write(r.content)
            print("   Download complete")
        else:
            print(f"Skipping {modality} model: already downloaded to {local_path}")


if __name__ == "__main__":
    main()
