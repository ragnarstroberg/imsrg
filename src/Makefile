CC = g++

#FLAGS = -larmadillo -O3 -std=c++11
FLAGS = -larmadillo -O3 -std=c++11  -lgsl #-fopenmp 
#FLAGS = -larmadillo -O3 -std=c++11 -pg

all: mysrg
	
OBJ = ModelSpace.o ReadWrite.o Operator.o HartreeFock.o IMSRGSolver.o AngMom.o imsrg_util.o

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(FLAGS) 

%.o: %.cc
	$(CC) -c $^ -o $@ $(FLAGS) 

clean:
	rm -f *.o
