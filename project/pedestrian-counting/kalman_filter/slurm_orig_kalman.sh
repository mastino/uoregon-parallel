#!/bin/bash
#
#SBATCH --job-name=kalman_trash
#SBATCH --output=data_orig_kalman/out_orig_kalman.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_orig_kalman/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_orig_kalman

tau_exec -ebs -T papi,openmp ./../orig_kalman_experiments.out 1000
