#!/bin/sh

module load orca/4.0.1
module load gcc/7.2.0
module unload openmpi/3.0.2
module unload openmpi/3.1.1
module unload openmpi/3.1.2
module load openmpi/2.0.4
/shared/apps/orca/orca_4_0_1_linux_x86-64_openmpi202/orca "$FILE_PATH.inp" > "$FILE_PATH.log" 