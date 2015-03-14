CC = g++

#FLAGS = -larmadillo -O3 -std=c++11
#FLAGS = -larmadillo -O3 -std=c++11  -lgsl -fPIC -fopenmp 
INCLUDE = -I./armadillo
#LIBS = -lblas -llapack -lgsl
LIBS = -lopenblas -llapack -lgsl
FLAGS =  -O3 -std=c++11  -fPIC -fopenmp # -DNO_ODE=1
#FLAGS = -larmadillo -O3 -std=c++11 -lgsl -fPIC # -v -da -Q

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

all: pyIMSRG.so libIMSRG.so
#all: mysrg pyIMSRG.so
	
#OBJ = ModelSpace.o ReadWrite.o Operator.o HartreeFock.o IMSRGSolver.o AngMom.o imsrg_util.o Operator3.o
OBJ = ModelSpace.o ReadWrite.o Operator.o HartreeFock.o IMSRGSolver.o AngMom.o imsrg_util.o

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

%.o: %.cc %.hh
	$(CC) -c $*.cc -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

libIMSRG.so: $(OBJ)
	$(CC) -shared -o $@ $^

python: pyIMSRG.so

pyIMSRG.so: $(OBJ) pyIMSRG.o
	$(CC) $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) -shared $(PYTHONFLAGS)

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c $^ -o $@ -fPIC $(FLAGS) $(INCLUDE) $(LIBS) $(PYTHONFLAGS)

clean:
	rm -f *.o *.so
