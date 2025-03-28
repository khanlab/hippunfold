import json

with open(snakemake.output.json, "w") as f:
    json.dump(snakemake.params.template_description, f, indent=4)
