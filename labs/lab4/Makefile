all: pass_serial pass_omp pass_tbb pass_cilk

pass_omp:
	

pass_tbb:
	
	
pass_cilk:
	

pass_serial: pass_serial.c
	icc -std=c99 -o pass_serial pass_serial.c -lcrypto

check_serial: pass_serial
	@./pass_serial 1120bfa29a3a834428dc477b68e2c872

check_omp: pass_omp
	

check_cilk: pass_cilk 
	

check_tbb: pass_tbb
	

clean:
	rm -f *.o pass_serial pass_omp pass_tbb pass_cilk
