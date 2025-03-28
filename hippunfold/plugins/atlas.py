from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import attrs
from git import GitCommandError, Repo
from snakebids import bidsapp
from snakebids.bidsapp.args import ArgumentGroups
from snakebids.plugins.base import PluginBase
from appdirs import AppDirs

logger = logging.getLogger(__name__)

import json
import os

# Global variable to store the commit hash
ATLAS_REPO_COMMIT = "7c292113c669a629773a644f82485ffde1842f52"


def get_download_dir():
    if "HIPPUNFOLD_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["HIPPUNFOLD_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("hippunfold", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir


def sync_atlas_repo():
    """
    Ensures the atlas folder is synced from the public GitHub repository using GitPython.
    """
    repo_url = "https://github.com/khanlab/hippunfold-atlases.git"
    atlas_dir = Path(get_download_dir()) / "hippunfold-atlases"

    try:
        if atlas_dir.exists() and (atlas_dir / ".git").exists():
            # If the directory exists and is a git repo, update it
            repo = Repo(atlas_dir)
            repo.remotes.origin.pull()
            repo.git.checkout(ATLAS_REPO_COMMIT)
        else:
            # If the directory does not exist, clone the repo
            repo = Repo.clone_from(repo_url, atlas_dir)
            repo.git.checkout(ATLAS_REPO_COMMIT)
    except GitCommandError as e:
        logger.info(f"Error syncing atlas repository: {e}")


def load_atlas_configs(atlas_dirs):
    """
    Loads surface atlas configurations from multiple directories.

    Args:
        atlas_dirs (list of str or Path): List of directories to search for surface atlas.
        The later directories in the list override the earlier ones.

    Returns:
        dict: A dictionary where keys are atlas names and values are their configurations.
    """
    atlas = {}

    for atlas_dir in atlas_dirs:
        atlas_path = Path(atlas_dir)
        if not atlas_path.exists():
            continue

        for subdir in atlas_path.iterdir():
            if subdir.is_dir():
                template_json = subdir / "template_description.json"
                if template_json.exists():
                    try:
                        with open(template_json, "r") as f:
                            config_data = json.load(f)
                        atlas[subdir.name.removeprefix("tpl-")] = (
                            config_data  # Override existing atlas if needed
                        )
                    except (json.JSONDecodeError, IOError) as e:
                        print(f"Warning: Failed to load {template_json}: {e}")

    return atlas


def get_atlas_configs():
    """
    Gets surface atlas configurations from both the centralized repo and cache directory.

    Returns:
        dict: Merged atlas configurations, prioritizing user-defined ones.
    """

    # Sync the atlas repository
    sync_atlas_repo()

    # Define search locations
    cache_dir = get_download_dir()

    atlas_dirs = []
    atlas_dirs.append(Path(cache_dir) / "hippunfold-atlases")
    #    atlas_dirs.append(Path(cache_dir) / "atlas_user")

    return load_atlas_configs(atlas_dirs)


@attrs.define
class AtlasConfig(PluginBase):
    """Dynamically add CLI parameters for the atlas.

    Parameters
    ----------
    argument_group
        Specify title of the group to which arguments should be added


    CLI Arguments
    ~~~~~~~~~~~~~
    Adds a atlas argument to the CLI, along with the related configuration
    for the atlas into the config.
    """

    argument_group: str | None = None
    atlas_config: dict = attrs.field(factory=get_atlas_configs)

    @bidsapp.hookimpl
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, argument_groups: ArgumentGroups
    ):
        """Add atlas parameters."""
        group = (
            argument_groups[self.argument_group]
            if self.argument_group is not None
            else parser
        )
        self.try_add_argument(
            group,
            "--atlas",
            choices=list(self.atlas_config.keys()),
            action="store",
            type=str,
            dest="atlas",
            default="multihist7",
            help=(
                "Select the atlas (unfolded space) to use for subfield labels. (default: %(default)s)"
            ),
        )

        self.try_add_argument(
            group,
            "--new_atlas_name",
            "--new-atlas-name",
            type=str,
            dest="new_atlas_name",
            help=("Name to use for the atlas created with group_create_atlas."),
        )
        self.try_add_argument(
            group,
            "--atlas-metrics",
            "--atlas_metrics",
            choices=["curvature", "gyrification", "thickness", "subfields"],
            action="store",
            type=str,
            dest="atlas_metrics",
            default=["curvature", "gyrification", "thickness"],
            nargs="+",
            help=(
                "Surface metrics to use when creating new atlas (default: %(default)s)"
            ),
        )
        self.try_add_argument(
            group,
            "--density",
            "--density",
            action="store",
            type=str,
            dest="density",
            default=[
                "hipp=272,1k,2k,4k,7k,10k,14k,18k,22k,27k,62k",
                "dentate=34,129,301,536,842,1k,2k,3k,8k,13k",
            ],
            nargs="+",
            help=(
                "Density for surface meshes. Use --atlas-list-density to get a listing of available densities (default: %(default)s)"
            ),
        )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Add atlas to config."""
        atlas = self.pop(namespace, "atlas")
        new_atlas_name = self.pop(namespace, "new_atlas_name")

        if (
            namespace["analysis_level"] == "group_create_atlas"
            and new_atlas_name == None
        ):
            raise argparse.ArgumentTypeError(
                "--new_atlas_name must be specified when using group_create_atlas"
            )

        config["atlas"] = atlas
        config["new_atlas_name"] = new_atlas_name
        config["atlas_metadata"] = self.atlas_config

        # parsing for density
        density_args = self.pop(namespace, "density")

        density = {}
        for d_arg in density_args:
            if "=" not in d_arg:
                raise argparse.ArgumentTypeError(
                    f"Invalid format: '{item}'. Expected format label=density1,density2,..."
                )

            label, values = d_arg.split("=", 1)
            if label not in ["hipp", "dentate"]:
                raise argparse.ArgumentTypeError(
                    f"Invalid label: '{label}'. Valid labels are hipp or dentate"
                )

            density[label] = values.split(",")

            # TODO: check if they are in the atlas

        config["density"] = density
