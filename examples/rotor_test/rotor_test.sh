#!/bin/sh
#set a job name
#SBATCH --job-name=rotor_test

#a file for job output, you can check job progress
#SBATCH --output=rotor_test.%a.slurm.log

# a file for errors from the job
#SBATCH --error=rotor_test.%a.slurm.log

#time you think you need; default is one day
# d-hh:mm:ss


#number of tasks you are requesting
#SBATCH -N 1
#SBATCH -n 10
##SBATCH --ntasks-per-node=2
##SBATCH --exclusive

#partition to use
#SBATCH --partition=ser-par-10g-2

#number of nodes to distribute n tasks across
#SBATCH -N 1

#an array job
#SBATCH --array=1-6


#####################################################

source activate rmg_env


echo $SLURM_ARRAY_TASK_ID
cd /gss_gpfs_scratch/harms.n/rotor_test # Edit this as need be
# the "stdbuf -o0 -e0"  and the "-u" are to disable buffering,
# so that you see output from the script in the log files immediately.
stdbuf -o0 -e0 python -u /home/harms.n/Code/AutoTST/examples/rotor_test/rotor_test.py > /home/harms.n/Code/AutoTST/examples/rotor_test/rotor_test.$SLURM_ARRAY_TASK_ID.combined.log 2>&1
