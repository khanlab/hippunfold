import json
import subprocess

template_description = dict(snakemake.params.template_description)

# record the ANTs version used for the template construction
try:
    ants_version = subprocess.run(
        ["antsRegistration", "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    ants_version = ants_version.removeprefix("ANTs Version:").strip()
except (OSError, subprocess.CalledProcessError):
    ants_version = "unknown"
template_description["GeneratedBy"] = template_description.get("GeneratedBy", []) + [
    {"Name": "ANTs", "Version": ants_version}
]

with open(snakemake.output.json, "w") as f:
    json.dump(template_description, f, indent=4)
