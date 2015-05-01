all: starhole_serial starhole_cilk starhole_openmp

starhole_serial: starhole_serial.cpp
	icpc -g -std=c++11 -o starhole_serial starhole_serial.cpp -Wall -Wextra -lopencv_core -lopencv_highgui -lm

starhole_cilk: starhole_cilk.cpp
	icpc -g -std=c++11 -o starhole_cilk starhole_cilk.cpp -Wall -Wextra -lopencv_core -lopencv_highgui -lm -lcilkrts

starhole_openmp: starhole_openmp.cpp
	icpc -g -std=c++11 -o starhole_openmp starhole_openmp.cpp -Wall -Wextra -lopencv_core -lopencv_highgui -lm -fopenmp

check_serial: starhole_serial
	time ./starhole_serial 10000 100  10 100 0 -100 0 0 100 0 -100
#	time ./starhole_serial 10000 100 200 100 0 -100 0 0 100 0 -100

clean:
	rm -f *.o starhole_serial starhole_openmp starhole_cilk out.jpg

.PHONY: clean
