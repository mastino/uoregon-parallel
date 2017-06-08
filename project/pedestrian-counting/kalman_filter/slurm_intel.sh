#!/bin/bash
#
#SBATCH --job-name=kalman_intel
#SBATCH --output=data_intel/out_intel.dat
#
#SBATCH --ntasks=1
#SBATCH -o data_intel/slurmjob-%j

module purge
module load intel/17
module load tau/intel
module load papi/5.5.1

export TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC

cd data_intel

tau_exec -ebs -T papi,openmp ./../run_intel_experiments.out 6 100000
