# 1120bfa29a3a834428dc477b68e2c872 -> 00000019
# 1120bfa29a3a834428dc477b68e77777 -> fail

CC=icc -std=c99


all: pass_serial pass_omp pass_tbb pass_cilk

pass_omp: pass_omp.c
	$(CC) -o pass_omp pass_omp.c -fopenmp  -lcrypto
	
	
pass_cilk: pass_cilk.c
	$(CC) -o pass_cilk pass_cilk.c -lcilkrts -lcrypto

pass_tbb: pass_tbb.cpp
	icpc -std=c++11 -fPIC -I$(TBBROOT)/include -o pass_tbb pass_tbb.cpp -L$(TBBLIB) -ltbb -lcrypto

pass_serial: pass_serial.c
	$(CC) -o pass_serial pass_serial.c -lcrypto

check_serial: pass_serial
<<<<<<< HEAD
	./pass_serial 1120bfa29a3a834428dc477b68e2c872
	./pass_serial 1120bfa29a3a834428dc477b68e77777
=======
	@./pass_serial 1120bfa29a3a834428dc477b68e2c872
>>>>>>> 8381863b2ac0b87cb59176985e2dd5bdfdafb1b8

check_omp: pass_omp
	./pass_omp 1120bfa29a3a834428dc477b68e2c872
	./pass_omp 1120bfa29a3a834428dc477b68e77777

check_cilk: pass_cilk 
	./pass_cilk 1120bfa29a3a834428dc477b68e2c872
	./pass_cilk 1120bfa29a3a834428dc477b68e77777
	

check_tbb: pass_tbb
	./pass_tbb 1120bfa29a3a834428dc477b68e2c872
	./pass_tbb 1120bfa29a3a834428dc477b68e77777
	

clean:
	rm -f *.o pass_serial pass_omp pass_tbb pass_cilk
