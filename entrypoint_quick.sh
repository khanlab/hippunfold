#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate snakebids-env
exec /src/hippunfold/run_quick.py "$@" 
