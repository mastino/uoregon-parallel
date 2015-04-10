CC=gcc
CC=tau_cc.sh -tau_makefile=/usr/local/packages/tau/x86_64/lib/Makefile.tau-icpc-papi-ompt-pdt-openmp
CFLAGS=-g -fopenmp -O3

all: tauTest

tauTest: tauTest.c
	$(CC) -o tauTest tauTest.c $(CFLAGS)

clean:
	rm -rf *.o tauTest  MULT* profile*
