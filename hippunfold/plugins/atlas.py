from __future__ import annotations

import argparse
import logging
import socket
import subprocess
from pathlib import Path
from typing import Any

import attrs
from appdirs import AppDirs
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

ATLASES = ["multihist7"]
DEFAULT_ATLAS = "multihist7"

ATLAS_METADATA = {
    "multihist7": {
        "metric_wildcards": ["curvature", "gyrification", "thickness"],
        "label_wildcards": ["hipp", "dentate"],
    }
}


DEFAULT_OUTPUT_DENSITY = ["8k"]
OUTPUT_DENSITIES = ["512", "2k", "8k", "18k"]

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


def get_unused_densities(output_density):
    """Gets the list of densities not used, so we can delete intermediate files."""
    return list(set(OUTPUT_DENSITIES) - set(output_density))


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
            choices=ATLASES,
            action="store",
            type=str,
            dest="atlas",
            default=DEFAULT_ATLAS,
            help=(
                "Select the atlas (unfolded space) to use for subfield labels. (default: %(default)s)"
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
            choices=OUTPUT_DENSITIES,
            nargs="+",
            help="Sets the output vertex density for participant-level results. Note: the density refers to the number of vertices in the hipp surface; the dentate has 1/4 the number of vertices.\n"
            + " (default: %(default)s)",
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
                "Method for defining subfields for the new atlas, either 'unfoldreg' to use unfoldreg with existing --atlas, or 'native' to use native space subfield segmentations. Note: if 'native' subfields are selected, data from both hemispheres will be concatenated."
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

        config["atlas"] = atlas
        config["new_atlas_name"] = new_atlas_name
        config["new_atlas_subfields_from"] = new_atlas_subfields_from
        config["output_density"] = output_density
        config["unfoldreg_density"] = OUTPUT_DENSITIES[-1]
        config["resample_factors"] = DEFAULT_RESAMPLE_FACTORS
        config["density_choices"] = resample_factors_to_densities(
            DEFAULT_RESAMPLE_FACTORS
        )
        config["unused_density"] = get_unused_densities(output_density)
        config["atlas_metadata"] = ATLAS_METADATA
