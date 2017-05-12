CC = icc -std=c99


all: crypto_serial crypto_omp crypto_tbb crypto_cilk

crypto_omp:
	$(CC) -o crypto_omp crypto_omp.c key.c -fopenmp

crypto_tbb:

crypto_cilk:

crypto_serial: serial.c key.c
	$(CC) -o crypto_serial serial.c key.c

check_serial: crypto_serial
	./crypto_serial superSecretText keys/key73 keys/key173 > cipher_serial &&\
	./crypto_serial cipher_serial   keys/key73 keys/key173

check_omp: crypto_omp
	./crypto_omp superSecretText keys/key73 keys/key173 > cipher_omp.encode &&\
	./crypto_omp cipher_omp.encode   keys/key73 keys/key173

check_tbb: crypto_tbb

check_cilk: crypto_cilk

clean:
	rm -f *.o crypto_serial crypto_omp crypto_tbb crypto_cilk cipher_serial
	
.PHONY:	clean
