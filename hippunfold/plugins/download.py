from __future__ import annotations

from typing import Any

import attrs
from snakebids import bidsapp
from snakebids.plugins.base import PluginBase

# last rules in each chain that populates the cache dir (HIPPUNFOLD_CACHE_DIR);
# their upstream download rules (download_surf_template_atlas,
# download_nnunet_model) are pulled in as dependencies. Only list the last rule of
# each chain: passing both a rule and its upstream rule to --until makes snakemake
# stop after the first job (snakemake/snakemake#823).
DOWNLOAD_RULES = [
    "download_extract_template",
    "cp_atlas_surf_gii",
    "cp_atlas_metric_gii",
]

# only defined when nnunet.smk is included (i.e. modality is not dsegtissue)
NNUNET_DOWNLOAD_RULES = [
    "unpack_nnunet_model",
]


@attrs.define
class DownloadConfig(PluginBase):
    """Restrict the ``download`` analysis level to the download rules.

    The ``download`` analysis level targets the same rule as ``participant``, but
    adds ``--until`` so only the jobs that populate the cache dir are run. Running
    this once before launching parallel participant-level runs (e.g. one per
    ``--participant-label``) avoids race conditions on the shared cache.
    """

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]):
        if config.get("analysis_level") != "download":
            return

        until = list(DOWNLOAD_RULES)
        if config.get("modality") != "dsegtissue":
            until.extend(NNUNET_DOWNLOAD_RULES)

        config["snakemake_args"] = [
            *config.get("snakemake_args", []),
            "--until",
            *until,
        ]
