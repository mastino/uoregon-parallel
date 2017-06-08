CC=taucc -std=c99 -tau:makefile $(TAU_DIR)/x86_64/lib/Makefile.tau-icpc-papi-ompt-pdt-openmp
CFLAGS= -g -fopenmp -fPIC
SFLAGS = --export=TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC
N = 1000
I = 10


test_linear_serial: linear_algebra.c
		$(CC) test_linear_algebra.c linear_algebra.c -o test_linear.out

test_linear_omp: linear_algebra_omp.c
		$(CC) test_linear_algebra.c linear_algebra_omp.c -o test_omp_linear.out $(CFLAGS)

test_linear_cilk: linear_algebra_cilk.c
		$(CC) test_linear_algebra.c -o test_linear.out linear_algebra_cilk.o 



linear_experiments: run_omp run_serial run_intel

serial_experiments: linear_algebra.c
		$(CC) run_linear_experiments.c linear_algebra.c -o run_serial_experiments.out $(CFLAGS) 

omp_experiments: linear_algebra_omp.c
		$(CC) run_linear_experiments.c linear_algebra_omp.c -o run_omp_experiments.out $(CFLAGS) 

intel_experiments: linear_algebra.c
		$(CC) -parallel run_linear_experiments.c linear_algebra.c -o run_intel_experiments.out $(CFLAGS) 

run_serial: serial_experiments
		sbatch slurm_serial.sh

run_omp: omp_experiments
		sbatch slurm_omp.sh

run_intel: intel_experiments
		sbatch slurm_intel.sh


orig_kalman_experiments: run_kalman_experiments.c
		$(CC) run_kalman_experiments.c kalman_filter.c linear_algebra.c -o orig_kalman_experiments.out $(CFLAGS) 

omp_kalman_experiments: run_kalman_experiments.c
		$(CC) run_kalman_experiments.c kalman_filter_omp.c linear_algebra_omp.c -o omp_kalman_experiments.out $(CFLAGS) 

serial_kalman_experiments: run_kalman_experiments.c
		$(CC) run_kalman_experiments.c kalman_filter_omp.c linear_algebra.c -o serial_kalman_experiments.out $(CFLAGS) 

intel_kalman_experiments: run_kalman_experiments.c
		$(CC) -parallel linear_algebra.c -c -o intel_linear.o $(CFLAGS) 
		$(CC) run_kalman_experiments.c kalman_filter_omp.c intel_linear.o -o intel_kalman_experiments.out $(CFLAGS) 

run_kalman: intel_kalman_experiments serial_kalman_experiments omp_kalman_experiments orig_kalman_experiments
		sbatch slurm_orig_kalman.sh
		sbatch slurm_serial_kalman.sh
		sbatch slurm_intel_kalman.sh
		sbatch slurm_omp_kalman.sh



test_kalman:
		$(CC) test-kalman.c kalman_filter.c linear_algebra.c -o test_kalman.out

test_kalman_omp:
		$(CC) test-kalman.c kalman_filter.c linear_algebra.c -o test_kalman_omp.out $(CFLAGS) 

test_kalman_db:
		$(CC) test-kalman.c kalman_filter.c linear_algebra.c -g -o test_kalman.out

linear_algebra_cilk.o: linear_algebra_cilk.c
	$(CC) -c linear_algebra_cilk.c -lcilkrts -lcrypto

linear_algebra_omp.o: linear_algebra_omp.c
	$(CC) -c linear_algebra_omp.c -fopenmp	

clean:
		rm -f *.o *.out slurmjob-* *.dat data_serial/slurmjob-* data_omp/slurmjob-* data_intel/slurmjob-*
		rm -f data_serial_kalman/slurmjob-* data_omp_kalman/slurmjob-* data_intel_kalman/slurmjob-* data_orig_kalman/slurmjob-*
		rm -f profile.*
		rm -rf MULTI__*
		rm -rf data_serial/MULTI__*
		rm -rf data_omp/MULTI__*
		rm -rf data_intel/MULTI__*
