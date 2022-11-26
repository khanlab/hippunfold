#!/bin/bash

set -euox pipefail

git restore ../README.md
tmpfile=`mktemp`
cp ../README.md $tmpfile
tail -n +5 $tmpfile > ../README.md
make latexpdf
cp _build/latex/hippunfold.pdf hippunfold_manual.pdf
git restore ../README.md

