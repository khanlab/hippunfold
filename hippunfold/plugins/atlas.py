from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import attrs
from appdirs import AppDirs
from git import GitCommandError, Repo
from snakebids import bidsapp
from snakebids.bidsapp.args import ArgumentGroups
from snakebids.plugins.base import PluginBase

logger = logging.getLogger(__name__)

import json
import os

try:
    from hippunfold.workflow.lib import utils as utils
except ImportError:
    from workflow.lib import utils as utils

# ====================================================================================
# This section is edited by hand:
# ------------------------------------------------------------------------------------
# Global variable to store the commit hash or branch name
ATLAS_REPO_COMMIT = "atlas-cli"

DEFAULT_OUTPUT_DENSITY = ["8k"]

# Default settings for atlas creation
DEFAULT_RESAMPLE_FACTORS = [
    12.5,
    25,
    50,
    75,
]  # percent, relative to native

# Default associated help (indicating approx vertex spacing for each factor)
DEFAULT_RESAMPLE_FACTORS_SPACING_HELP = ", ".join(["~2mm", "~1mm", "~0.5mm", "0.3mm"])
# ====================================================================================


# helper functions for resample factors and density


def resample_to_density(r):
    nverts = 256 * 128 * (r / 100) ** 2
    if nverts > 800:
        return str(round(nverts / 1000)) + "k"
    else:
        return str(int(nverts))


def resample_factors_to_densities(resample_factors):

    return [resample_to_density(r) for r in resample_factors]


def sort_densities(densities):
    def parse_density(d):
        if isinstance(d, str) and d.endswith("k"):
            return int(float(d[:-1]) * 1000)
        return int(d)

    return sorted(densities, key=parse_density)


def get_all_densities(atlas_config):
    """Return a sorted list of all unique output densities across all atlases."""
    densities = set()
    for config in atlas_config.values():
        densities.update(config.get("density_wildcards", []))
    return ["native"] + sort_densities(densities)


def get_unfoldreg_density(atlas_config, atlas):
    """Get the density to use for unfoldreg, which is the highest density available for the chosen atlas."""
    return sort_densities(atlas_config[atlas]["density_wildcards"])[-1]


def get_unused_densities(atlas_config, atlas, output_density):
    """Gets the list of densities not used, so we can delete intermediate files."""
    return list(
        set(["native"] + atlas_config[atlas]["density_wildcards"]) - set(output_density)
    )


def format_density_help(atlas_config):
    lines = ["Available output densities per atlas:\n"]
    for atlas, config in atlas_config.items():
        if "density_wildcards" in config:
            densities = ", ".join(config["density_wildcards"])
            lines.append(f"  {atlas}=[{densities}]")
    return "\n".join(lines)


def validate_output_density(atlas, output_densities, atlas_config):
    """
    Validate that each output_density is allowed for the selected atlas.

    Parameters:
        atlas (str): The name of the selected atlas.
        output_densities (str or list): One or more densities to validate.
        atlas_config (dict): Dictionary of atlas options, with allowed densities under ['density_wildcards'].

    Raises:
        ValueError: If any of the provided densities are invalid.
    """
    if isinstance(output_densities, str):
        output_densities = [output_densities]

    allowed = set(atlas_config[atlas]["density_wildcards"]) | {"native"}

    invalid = [d for d in output_densities if d not in allowed]
    if invalid:
        raise ValueError(
            f"Invalid output_density value(s) for atlas '{atlas}': {invalid}. "
            f"Allowed values: {sorted(allowed)}"
        )


# helper functions for hippunfold-atlases
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

    return load_atlas_configs(atlas_dirs)


# snakebids plugin definition


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
            "--new_atlas_subfields_from",
            "--new-atlas-subfields-from",
            type=str,
            choices=["unfoldreg", "native"],
            dest="new_atlas_subfields_from",
            help=(
                "Method for defining subfields for the new atlas, either 'unfoldreg' to use unfoldreg with existing --atlas, or 'native' to use native space subfield segmentations."
            ),
        )
        self.try_add_argument(
            group,
            "--new-atlas-metrics",
            "--new_atlas_metrics",
            choices=["curvature", "gyrification", "thickness", "myelin"],
            action="store",
            type=str,
            dest="new_atlas_metrics",
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
            default=DEFAULT_OUTPUT_DENSITY,
            choices=get_all_densities(self.atlas_config),
            nargs="+",
            help="Sets the output vertex density for participant-level results. Note: the density refers to the number of vertices in the hipp surface; the dentate has 1/4 the number of vertices.\n"
            + format_density_help(self.atlas_config)
            + " (default: %(default)s)",
        )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Add atlas to config."""
        atlas = self.pop(namespace, "atlas")
        new_atlas_name = self.pop(namespace, "new_atlas_name")
        new_atlas_subfields_from = self.pop(namespace, "new_atlas_subfields_from")
        output_density = self.pop(namespace, "output_density")

        if namespace["analysis_level"] == "group_create_atlas":

            if new_atlas_name == None:
                raise argparse.ArgumentTypeError(
                    "--new_atlas_name must be specified when using group_create_atlas"
                )
            if new_atlas_subfields_from == None:
                raise argparse.ArgumentTypeError(
                    "--new_atlas_subfields_from must be specified when using group_create_atlas"
                )

        validate_output_density(atlas, output_density, self.atlas_config)

        config["atlas"] = atlas
        config["new_atlas_name"] = new_atlas_name
        config["new_atlas_subfields_from"] = new_atlas_subfields_from
        config["atlas_metadata"] = self.atlas_config
        config["output_density"] = output_density
        config["unfoldreg_density"] = get_unfoldreg_density(self.atlas_config, atlas)
        config["resample_factors"] = DEFAULT_RESAMPLE_FACTORS
        config["density_choices"] = resample_factors_to_densities(
            DEFAULT_RESAMPLE_FACTORS
        )
        config["unused_density"] = get_unused_densities(
            self.atlas_config, atlas, output_density
        )
