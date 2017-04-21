CC=icc
CXX=icpc
CFLAGS=-g -O3
CXXFLAGS=-g -O3

all: openmpNotQuiteRight cilkNotQuiteRight tbbNotQuiteRight


openmpNotQuiteRight: openmpNotQuiteRight.c
	$(CC) $(CFLAGS) -o openmpNotQuiteRight openmpNotQuiteRight.c -fopenmp

cilkNotQuiteRight: cilkNotQuiteRight.c
	$(CC) $(CFLAGS) -o cilkNotQuiteRight cilkNotQuiteRight.c -lcilkrts

tbbNotQuiteRight: tbbNotQuiteRight.cpp
	$(CXX) $(CXXFLAGS) -o tbbNotQuiteRight tbbNotQuiteRight.cpp -ltbb

check: cilkNotQuiteRight
	@cilkNotQuiteRight 4

clean: 
	rm tbbNotQuiteRight cilkNotQuiteRight openmpNotQuiteRight
