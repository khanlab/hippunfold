#!/usr/bin/env python3
import os
from pathlib import Path

from snakebids import bidsapp, plugins

try:
    from hippunfold.plugins import atlas as atlas_plugin  # Works when run as a package
    from hippunfold.workflow.lib import utils as utils
except ImportError:
    from plugins import atlas as atlas_plugin  # Works when run directly
    from workflow.lib import utils as utils


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


# Set the environment variable SNAKEMAKE_CONDA_PREFIX if not already set
if not "SNAKEMAKE_CONDA_PREFIX" in os.environ:
    os.environ["SNAKEMAKE_CONDA_PREFIX"] = str(Path(utils.get_download_dir()) / "conda")


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
