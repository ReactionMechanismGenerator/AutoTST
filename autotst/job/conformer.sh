#!/bin/sh
conda activate $ENV || source activate $ENV
python $AUTOTST/autotst/job/conformer.py "$SMILES" "$DIRECTORY"