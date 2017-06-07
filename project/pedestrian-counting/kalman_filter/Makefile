CC=taucc -std=c99 -tau:makefile $(TAU_DIR)/x86_64/lib/Makefile.tau-icpc-papi-ompt-pdt-openmp
CFLAGS= -g -fopenmp -fPIC
SFLAGS = --export=TAU_METRICS=TIME,PAPI_TOT_INS,PAPI_TOT_CYC
N = 1000
I = 10

all: run_omp run_serial

test_linear_serial: linear_algebra.c
		$(CC) test_linear_algebra.c linear_algebra.c -o test_linear.out

test_linear_omp: linear_algebra_omp.c
		$(CC) test_linear_algebra.c linear_algebra_omp.c -o test_omp_linear.out $(CFLAGS)

test_linear_cilk: linear_algebra_cilk.c
		$(CC) test_linear_algebra.c -o test_linear.out linear_algebra_cilk.o 


serial_experiments: linear_algebra.c
		$(CC) run_linear_experiments.c linear_algebra.c -o run_serial_experiments.out $(CFLAGS) 

omp_experiments: linear_algebra_omp.c
		$(CC) run_linear_experiments.c linear_algebra_omp.c -o run_omp_experiments.out $(CFLAGS) 

cilk_experiments: linear_algebra_cilk.c
		$(CC) run_linear_experiments.c linear_algebra_cilk.c -o run_cilk_experiments.out $(CFLAGS) 


run_serial: serial_experiments
		sbatch slurm_serial.sh

run_omp: omp_experiments
		sbatch slurm_omp.sh

run_cilk: cilk_experiments
		sbatch slurm_cilk.sh


test_kalman:
		$(CC) test-kalman.c kalman_filter.c linear_algebra.c -o test_kalman.out

test_kalman_db:
		$(CC) test-kalman.c kalman_filter.c linear_algebra.c -g -o test_kalman.out

linear_algebra_cilk.o: linear_algebra_cilk.c
	$(CC) -c linear_algebra_cilk.c -lcilkrts -lcrypto

linear_algebra_omp.o: linear_algebra_omp.c
	$(CC) -c linear_algebra_omp.c -fopenmp	

clean:
		rm -f *.o *.out slurmjob-* *.dat data_serial/slurmjob-* data_omp/slurmjob-*
		rm -rf MULTI__*
		rm -rf data_serial/MULTI__*
		rm -rf data_omp/MULTI__*
