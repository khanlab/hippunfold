#!/bin/bash

# Run the command, passing the default flags and any additional arguments
exec conda run --no-capture-output -n snakebids-env /src/hippunfold/run.py "$@"
