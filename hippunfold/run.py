#!/usr/bin/env python3
import os
from pathlib import Path

from snakebids import bidsapp, plugins

from hippunfold.workflow.lib import utils as utils

try:
    from hippunfold.plugins import atlas as atlas_plugin  # Works when run as a package
except ImportError:
    from plugins import atlas as atlas_plugin  # Works when run directly


if "__file__" not in globals():
    __file__ = "../hippunfold/run.py"


app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="hippunfold"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
        atlas_plugin.AtlasConfig(argument_group="ATLASES"),
    ]
)

# Set the conda prefix directory
conda_prefix = str(utils.get_download_dir()) + "/" + "conda"

# Set the environment variable SNAKEMAKE_CONDA_PREFIX
os.environ["SNAKEMAKE_CONDA_PREFIX"] = str(conda_prefix)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
