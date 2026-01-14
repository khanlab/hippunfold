
import os
from snakebids.app import SnakeBidsApp
from snakebids.cli import add_dynamic_args
import json

# Mimic the app structure using the existing config
app = SnakeBidsApp(
    os.path.abspath("hippunfold"),
    configfile_path="config/snakebids.yml",
)

# Parse args simluating the user command
# hippunfold test_data/bids_dhcp/ test_output/bids_dhcp participant --modality T2w --force_nnunet_model neonate_synthseg --template dHCP --cores 8
args = [
    "test_data/bids_dhcp/",
    "test_output/bids_dhcp",
    "participant",
    "--modality", "T2w",
    "--force_nnunet_model", "neonate_synthseg",
    "--force_nnunet_model", "neonate_synthseg",
    "--template", "dHCP",
]

# Simulate what happens inside run_snakemake (or rather before it, since we want to check config)
# Actually, the parsing happens when accessing `app.parser`.

# add_dynamic_args(app.parser, app.config["parse_args"], app.config["pybids_inputs"])
parsed_args = app.parser.parse_args(args)

# SnakeBidsApp.update_config update config with parsed args
app.config.update(vars(parsed_args))

# Now the critical part: generate_inputs
# internal method usually called by SnakeBids
from snakebids.core.input_generation import generate_inputs

# config should have bids_dir and pybids_inputs
# snakebids.yml has T2w wildcards: subject, session, acquisition, run
# Attempt to simplify them to match the file
print("Original T2w wildcards:", app.config["pybids_inputs"]["T2w"]["wildcards"])
app.config["pybids_inputs"]["T2w"]["wildcards"] = ["subject", "session"]
print("Modified T2w wildcards:", app.config["pybids_inputs"]["T2w"]["wildcards"])

print(f"BIDS DIR: {app.config['bids_dir']}")

inputs = generate_inputs(
    bids_dir=app.config["bids_dir"],
    pybids_inputs=app.config["pybids_inputs"],
    pybids_database_dir=app.config.get("pybids_database_dir"),
    pybids_reset_database=app.config.get("pybids_reset_database"),
    derivatives=app.config.get("derivatives"),
    participant_label=app.config.get("participant_label"),
    exclude_participant_label=app.config.get("exclude_participant_label"),
)

print("Inputs keys:", inputs.keys())
if "T2w" in inputs:
    print("T2w input files:", inputs["T2w"])
else:
    print("T2w NOT FOUND in inputs")

print("-" * 20)
print("Debugging with BIDSLayout directly")
from bids import BIDSLayout
layout = BIDSLayout(app.config['bids_dir'], validate=False) # validate=False often helps with non-standard datasets
print(f"Layout files count: {len(layout.get())}")
print("T2w files in layout (basic):", layout.get(suffix='T2w', extension='.nii.gz'))
print("T2w files in layout (space=None):", layout.get(suffix='T2w', extension='.nii.gz', space=None))
print("T2w files in layout (datatype=anat):", layout.get(suffix='T2w', extension='.nii.gz', datatype='anat'))

print("-" * 20)
print("Testing removing space: null from config")
if "space" in app.config["pybids_inputs"]["T2w"]["filters"]:
    del app.config["pybids_inputs"]["T2w"]["filters"]["space"]
    
inputs_nospace = generate_inputs(
    bids_dir=app.config["bids_dir"],
    pybids_inputs=app.config["pybids_inputs"],
)
if "T2w" in inputs_nospace:
    print("T2w found WITHOUT space filter:", inputs_nospace["T2w"])
else:
    print("T2w NOT found without space filter")

