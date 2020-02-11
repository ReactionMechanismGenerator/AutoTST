# Commands on Discovery
`singularity shell --bind /bin/ --bind /etc/slurm/ --bind /etc/munge --bind /var/log/munge/ --bind /usr/local/bin/ --bind /usr/local/sbin/ --bind /etc/profile.d/  --bind /usr/lib64/ --bind /etc/passwd  --bind /var/run/munge/ --bind /usr/share/Modules/  --bind /shared/centos7/ --home $PWD --writable-tmpfs  --bind /scratch/westgroup/singularity_images/autotst/  autotst_latest.sif`
`source /apps/.bashrc` 
`unset DFTB_PREFIX` 
`export DFTB_PREFIX=/apps/dftbplus-19.1.x86_64-linux/halorg-0-1/` 
`export AUTOTST=/apps/AutoTST/` 
`unset DFTB_COMMAND` 
`export DFTB_COMMAND=/apps/dftbplus-19.1.x86_64-linux/bin/dftb+` 
`source /shared/centos7/gaussian/g16/bsd/g16.profile` 
`export AUTOTST_SCR_DIR=/scratch/westgroup/singularity_images/autotst` 
`export PYTHONPATH=$AUTOTST:$PYTHONPATH`