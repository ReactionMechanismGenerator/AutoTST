#!/bin/sh
export GAUSS_SCRDIR="/scratch/$USER/$SLURM_JOBID"
$COMMAND "$FILE_PATH.com" > "$FILE_PATH.log"
