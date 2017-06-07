#!/bin/bash
#
#SBATCH --job-name=kalman_omp
#SBATCH --output=data_omp/out_omp.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_omp/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_omp

tau_exec -ebs -T papi,openmp ./../run_omp_experiments.out 6 100000
