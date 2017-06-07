#!/bin/bash
#
#SBATCH --job-name=kalman_cilk
#SBATCH --output=data_cilk/out_cilk.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_cilk/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_cilk

tau_exec -ebs -T papi,cilk,openmp ./../run_cilk_experiments.out 6 100000
