# 1120bfa29a3a834428dc477b68e2c872 -> 00000019
# 1120bfa29a3a834428dc477b68e77777 -> fail

CC=icc -std=c99


all: pass_serial pass_omp pass_tbb pass_cilk

pass_omp: pass_omp.c
	$(CC) -o pass_omp pass_omp.c -fopenmp  -lcrypto
	
	
pass_cilk: pass_cilk.c
	

pass_tbb: pass_tbb.c
	

pass_serial: pass_serial.c
	$(CC) -o pass_serial pass_serial.c -lcrypto

check_serial: pass_serial
	@pass_serial 1120bfa29a3a834428dc477b68e2c872

check_omp: pass_omp
	time ./pass_omp 1120bfa29a3a834428dc477b68e2c872
	time ./pass_omp 1120bfa29a3a834428dc477b68e77777
	time ./pass_serial 1120bfa29a3a834428dc477b68e2c872
	time ./pass_serial 1120bfa29a3a834428dc477b68e77777

check_cilk: pass_cilk 
	

check_tbb: pass_tbb
	

clean:
	rm -f *.o pass_serial pass_omp pass_tbb pass_cilk
