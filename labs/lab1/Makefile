# For TBB, make sure you first load the appropriate module to get 
# the values for TBBROOT and TBBLIB, e.g., module load intel/tbb
all: mpWorld cilkWorld tbbWorld

mpWorld: helloWorld.c
	icc -o mpWorld helloWorld.c -fopenmp

check:
	export OMP_NUM_THREADS=1 &&	time ./mpWorld
	export OMP_NUM_THREADS=2 &&	time ./mpWorld
	export OMP_NUM_THREADS=4 &&	time ./mpWorld
	export OMP_NUM_THREADS=8 &&	time ./mpWorld
	export OMP_NUM_THREADS=16 &&	time ./mpWorld
	export OMP_NUM_THREADS=32 &&	time ./mpWorld

cilkWorld: cilkWorld.c
	icc -o cilkWorld cilkWorld.c -lcilkrts

tbbWorld: tbbWorld.cpp
	icpc -fPIC -I$(TBBROOT)/include -o tbbWorld tbbWorld.cpp -L$(TBBLIB) -ltbb

clean: 
	rm tbbWorld cilkWorld mpWorld 
