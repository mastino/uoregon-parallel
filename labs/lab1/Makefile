# For TBB, make sure you first load the appropriate module to get 
# the values for TBBROOT and TBBLIB, e.g., module load intel/tbb
all: mpWorld cilkWorld tbbWorld

mpWorld: helloWorld.c
	icc -o mpWorld helloWorld.c -fopenmp

cilkWorld: cilkWorld.c
	icc -o cilkWorld cilkWorld.c -lcilkrts

tbbWorld: tbbWorld.cpp
	icpc -fPIC -I$(TBBROOT)/include -o tbbWorld tbbWorld.cpp -L$(TBBLIB) -ltbb

clean: 
	rm tbbWorld cilkWorld mpWorld 
