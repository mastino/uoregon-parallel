CC = icc -std=c99
KEYS = keys/key1 keys/key14 keys/key181 keys/key222 keys/key264 keys/key305 keys/key347 keys/key389 keys/key43 keys/key471 keys/key512 keys/key554 keys/key596 keys/key637 keys/key679 keys/key72 keys/key761 keys/key802 keys/key844 keys/key886 keys/key927 keys/key969

all: crypto_serial crypto_omp crypto_tbb crypto_cilk

crypto_omp: crypto_omp.c
	$(CC) -o crypto_omp crypto_omp.c key.c -fopenmp

crypto_tbb:
	icpc -std=c++11 -fPIC -I$(TBBROOT)/include -o crypto_tbb crypto_tbb.cpp key.c -L$(TBBLIB) -ltbb -lcrypto

crypto_cilk: crypto_cilk.c
	$(CC) -o crypto_cilk crypto_cilk.c key.c -lcilkrts -lcrypto

crypto_serial: serial.c key.c
	$(CC) -o crypto_serial serial.c key.c

check_serial: crypto_serial
	./crypto_serial 19033.txt $(KEYS) > cipher_serial &&\
	./crypto_serial cipher_serial $(KEYS)

check_omp: crypto_omp
	./crypto_omp 19033.txt $(KEYS)  > cipher_omp.encode &&\
	./crypto_omp cipher_omp.encode   $(KEYS) 

check_tbb: crypto_tbb
	./crypto_tbb 19033.txt $(KEYS) > cipher_tbb.encode &&\
	./crypto_tbb cipher_tbb.encode   $(KEYS)

check_cilk: crypto_cilk
	./crypto_cilk 19033.txt $(KEYS) > cipher_cilk.encode &&\
	./crypto_cilk cipher_cilk.encode   $(KEYS)

clean:
	rm -f *.o crypto_serial crypto_omp crypto_tbb crypto_cilk cipher_serial
	
.PHONY:	clean
