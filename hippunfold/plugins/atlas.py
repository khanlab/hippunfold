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


try:
    from hippunfold.workflow.lib import utils as utils
except ImportError:
    from workflow.lib import utils as utils


# Global variable to store the commit hash
ATLAS_REPO_COMMIT = "679f5d1525a82dbbd4327c265a15b5729a32f263"

RESAMPLING_FACTORS = [
    25,
    50,
    75,
]  # percent, relative to native
# in multihist7, this corresponds to ~1mm, ~0.5mm, ~0.33m vertex distances

# naming for native similar to fsnative
ATLAS_DENSITY_CHOICES = [
    "native",
]
# naming convention based on fsLR32k, etc.
# NOTE: this is based on the hipp surface not the dentate surface.
for R in RESAMPLING_FACTORS:
    ATLAS_DENSITY_CHOICES.append(str(int((256 * 128 * (R / 100) ** 2) / 1000)) + "k")
# See output file tpl-ATLAS_desc-resample2density_mapping.csv for estimates of vertex spacing in mm

ATLAS_DENSITY_DEFAULT = ATLAS_DENSITY_CHOICES[
    1
]  # also density that is used for unfoldreg, cannot set this to native


def sync_atlas_repo():
    """
    Ensures the atlas folder is synced from the public GitHub repository using GitPython.
    """
    repo_url = "https://github.com/khanlab/hippunfold-atlases.git"
    atlas_dir = Path(utils.get_download_dir()) / "hippunfold-atlases"

    try:
        if atlas_dir.exists() and (atlas_dir / ".git").exists():
            repo = Repo(atlas_dir)
            repo.git.fetch()  # Make sure latest commits are available
        else:
            repo = Repo.clone_from(repo_url, atlas_dir)

        # Explicitly checkout desired commit
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
    cache_dir = utils.get_download_dir()

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
            "--output-density",
            "--output_density",
            action="store",
            type=str,
            dest="output_density",
            default=[ATLAS_DENSITY_DEFAULT],
            choices=ATLAS_DENSITY_CHOICES,
            nargs="+",
            help=(
                "Sets the output vertex density for results, using the same vertex density for hipp and dentate (default: %(default)s)"
            ),
        )
        self.try_add_argument(
            group,
            "--resample-factors",
            "--resample_factors",
            action="store",
            type=str,
            dest="output_density",
            default=[RESAMPLING_FACTORS],
            choices=RESAMPLING_FACTORS,
            nargs="+",
            help=(
                "Sets the downsampling factors of the surface mesh relative to native (default: %(default)s)"
            ),
        )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Add atlas to config."""
        atlas = self.pop(namespace, "atlas")
        new_atlas_name = self.pop(namespace, "new_atlas_name")
        output_density = self.pop(namespace, "output_density")

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
        config["output_density"] = output_density
        config["unfoldreg_density"] = ATLAS_DENSITY_DEFAULT
        config["resample_factors"] = RESAMPLING_FACTORS
        config["density_choices"] = ATLAS_DENSITY_CHOICES
        config["unused_density"] = list(
            set(ATLAS_DENSITY_CHOICES) - set(output_density)
        )
