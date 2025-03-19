#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins

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


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
