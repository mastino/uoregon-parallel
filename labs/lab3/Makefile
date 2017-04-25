
CC=icc

# Instrument the code with TAU
CC=tau_cc.sh -tau_makefile=$(TAU_DIR)/x86_64/lib/Makefile.tau-icpc-papi-ompt-pdt-openmp
CFLAGS=-g -fopenmp -O3

all: tauTest

tauTest: tauTest.c
	$(CC) -o tauTest tauTest.c $(CFLAGS)

clean:
	rm -rf *.o tauTest  MULT* profile*
