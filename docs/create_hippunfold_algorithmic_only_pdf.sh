#!/bin/bash

set -euox pipefail

git restore index.md conf.py
cp  supp_methods_index.md index.md
cp  supp_methods_conf.py conf.py
make clean
make latexpdf
cp _build/latex/hippunfoldalgorithmicdetails.pdf hippunfold_algorithmic_details.pdf
git restore index.md conf.py

#git add hippunfold_manual.pdf
#git commit -m "updated hippunfold manual"

