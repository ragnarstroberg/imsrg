CC = g++
# Need to do a better job with the dependencies...

#FLAGS = -larmadillo -O3 -std=c++11
FLAGS = -larmadillo -O3 -std=c++11  -lgsl -fPIC -fopenmp 
#FLAGS = -larmadillo -O3 -std=c++11 -lgsl -fPIC # -v -da -Q

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

all: mysrg pyIMSRG.so
	
OBJ = ModelSpace.o ReadWrite.o Operator.o HartreeFock.o IMSRGSolver.o AngMom.o imsrg_util.o

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(FLAGS) 

%.o: %.cc
	$(CC) -c $^ -o $@ $(FLAGS) 



python: pyIMSRG.so

pyIMSRG.so: pyIMSRG.o $(OBJ)
	$(CC) $^ -o $@ $(FLAGS) -shared $(PYTHONFLAGS)

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c $^ -o $@ -fPIC $(FLAGS) $(PYTHONFLAGS)

clean:
	rm -f *.o
