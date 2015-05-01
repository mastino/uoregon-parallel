


cc=icpc
ccflags= -MMD -g -Wall
ldflags= -lopencv_core -lopencv_highgui -lopencv_imgproc

all: pipeline_serial

pipeline_serial.o: pipeline_serial.cpp
	$(cc) $(ccflags) -c pipeline_serial.cpp

pipeline_serial: pipeline_serial.o
	$(cc) $(ccflags) -o pipeline_serial pipeline_serial.o $(ldflags)

check_serial: pipeline_serial
	./pipeline_serial /tmp/vframe/*

pipeline_parallel:

check_parallel:


.PHONY: clean

clean:
	rm -f  *.o *~  pipeline_serial


