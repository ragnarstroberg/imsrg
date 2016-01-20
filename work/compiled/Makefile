
CC = g++

FLAGS = -std=c++11 -fopenmp  -O3
INCLUDE = -I$(HOME)/include -I$(HOME)/include/armadillo
#LIBS = -lIMSRG -lopenblas -llapack -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS = -L$(HOME)/lib -lIMSRG -lopenblas  -lgsl -lz -lhdf5 -lhdf5_cpp
#LIBS = -L$(HOME)/lib64 -L$(HOME)/lib -lopenblas -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS += -llapack
ifeq ($(DEBUG),on)
 FLAGS += -g
endif

all: exe
exe: SingleRef ValenceSpace

%:	%.cc
	$(CC) $(INCLUDE) $^ -o $@ $(FLAGS) $(LIBS)

