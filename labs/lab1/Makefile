
all: mpWorld cilkWorld tbbWorld

mpWorld: helloWorld.c
	icc -o mpWorld helloWorld.c -fopenmp

cilkWorld: cilkWorld.c
	icc -o cilkWorld cilkWorld.c -lcilkrts

tbbWorld: tbbWorld.cpp
	icpc -o tbbWorld tbbWorld.cpp -ltbb

clean: 
	rm tbbWorld cilkWorld mpWorld 
