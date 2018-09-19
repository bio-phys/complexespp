#!/bin/bash

# join stdout and stderr
#$ -j y
# change to current dir
#$ -cwd
# 24h Queue
#$ -l h_rt=24:00:00
#$ -N twoD-1q0w_A_B
# number of cores
#$ -pe impi_hydra 24


export OMP_PROC_BIND=TRUE
export OMP_WAIT_POLICY=PASSIVE

set -xe

complexespp=complexes++
$complexespp -c config --multidir=$(ls | grep T_) --replex=100 --replex-accept=remc --nb-threads=10 --backup=False 2>/dev/null
