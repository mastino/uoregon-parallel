test_linear:
		gcc test_linear_algebra.c linear_algebra.c -o test_linear.out

test_linear_omp: linear_algebra_omp.o
		gcc test_linear_algebra.c linear_algebra_omp.o -o test_omp_linear.out -fopenmp

run_omp_experiments:
		gcc run_linear_experiments.c linear_algebra_omp.o -o run_omp_experiments.out -fopenmp 

test_linear_cilk: linear_algebra_cilk.o
		gcc test_linear_algebra.c -o test_linear.out linear_algebra_cilk.o 

test_kalman:
		gcc test-kalman.c kalman_filter.c linear_algebra.c -o test_kalman.out

test_kalman_db:
		gcc test-kalman.c kalman_filter.c linear_algebra.c -g -o test_kalman.out

linear_algebra_cilk.o: linear_algebra_cilk.c
	$(CC) -c linear_algebra_cilk.c -lcilkrts -lcrypto

linear_algebra_omp.o: linear_algebra_omp.c
	$(CC) -c linear_algebra_omp.c -fopenmp	

clean:
		rm -f *.o *.out