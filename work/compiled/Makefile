
CC = g++
CCFLAGS = -std=c++11 -fopenmp  #-g
INCLUDE = -I$(HOME)/include -I$(HOME)/include/armadillo
LIBS = -L$(HOME)/lib64 -lIMSRG -lopenblas -llapack -lgsl -lz -lhdf5_cpp

all: exe
exe: SingleRef ValenceSpace

%:	%.cc
	$(CC) $(INCLUDE) $^ -o $@ $(CCFLAGS) $(LIBS)

