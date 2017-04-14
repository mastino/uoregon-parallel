CC=icc
CXX=icpc

all: openmpNotQuiteRight cilkNotQuiteRight tbbNotQuiteRight


openmpNotQuiteRight: openmpNotQuiteRight.c
	$(CC) -o openmpNotQuiteRight openmpNotQuiteRight.c -fopenmp

cilkNotQuiteRight: cilkNotQuiteRight.c
	$(CC) -o cilkNotQuiteRight cilkNotQuiteRight.c -lcilkrts

tbbNotQuiteRight: tbbNotQuiteRight.cpp
	$(CXX) -o tbbNotQuiteRight tbbNotQuiteRight.cpp -ltbb

check: cilkNotQuiteRight
	@cilkNotQuiteRight 4

clean: 
	rm tbbNotQuiteRight cilkNotQuiteRight openmpNotQuiteRight
