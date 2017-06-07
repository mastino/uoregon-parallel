#!/bin/bash
#
#SBATCH --job-name=kalman_serial
#SBATCH --output=data_serial/out_serial.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_serial/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_serial

tau_exec -ebs -T papi,openmp ./../run_serial_experiments.out 6 100000
