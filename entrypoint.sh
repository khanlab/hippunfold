#!/bin/bash

# Default flags to pass to run.py
DEFAULT_FLAGS="--use-conda --conda-prefix /src/conda-envs"

# Run the command, passing the default flags and any additional arguments
exec conda run --no-capture-output -n snakebids-env /src/hippunfold/run.py "$DEFAULT_FLAGS" "$@"
