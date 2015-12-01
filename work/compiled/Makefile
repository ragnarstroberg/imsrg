
CC = g++
CCFLAGS = -std=c++11 -fopenmp 
INCLUDE = -I$(HOME)/include -I$(HOME)/include/armadillo
LIBS = -L$(HOME)/lib64 -L$(HOME)/lib -lIMSRG -lopenblas -lgsl -lz -lhdf5 -lhdf5_cpp
#LIBS = -L$(HOME)/lib64 -L$(HOME)/lib -lopenblas -lgsl -lz -lhdf5 -lhdf5_cpp

all: exe
exe: SingleRef ValenceSpace

%:	%.cc
	$(CC) $(INCLUDE) $^ -o $@ $(CCFLAGS) $(LIBS)

