#!/bin/bash
#
#SBATCH --job-name=kalman_trash
#SBATCH --output=data_serial_kalman/out_serial_kalman.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_serial_kalman/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_serial_kalman

tau_exec -ebs -T papi,openmp ./../serial_kalman_experiments.out 1000
