#!/bin/sh
conda activate $ENV
python $AUTOTST/autotst/job/conformer.py "$SMILES" "$DIRECTORY"
