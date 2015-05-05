CC = g++

#FLAGS = -larmadillo -O3 -std=c++11
#FLAGS = -larmadillo -O3 -std=c++11  -lgsl -fPIC -fopenmp 
INCLUDE = -I./armadillo
#LIBS = -lblas -llapack -lgsl
LIBS = -lopenblas -llapack -lgsl
#FLAGS =  -O3 -std=c++11  -fPIC -fopenmp # -DNO_ODE=1
#FLAGS =  -std=c++11 -fopenmp -fPIC -g # -DNO_ODE=1
FLAGS =  -O2 -march=native -std=c++11 -fopenmp -fPIC -DNO_ODE=1 -g
#FLAGS =  -std=c++11 -fopenmp -fPIC -g  -DNO_ODE=1
SOFLABS = -std=c++11 -fopenmp -fPIC -g
#FLAGS = -larmadillo -O3 -std=c++11 -lgsl -fPIC # -v -da -Q
WARNFLAGS = -Wall -Wno-comment
FLAGS += $(WARNFLAGS)
SOFLAGS += $(WARNFLAGS)

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

all: pyIMSRG.so libIMSRG.so
#all: libIMSRG.so
#all: mysrg pyIMSRG.so
	
OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
      HartreeFock.o Generator.o IMSRGSolver.o AngMom.o imsrg_util.o
#OBJ =  CrapOperator.o 
      

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

%.o: %.cc %.hh
	$(CC) -c $*.cc -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

libIMSRG.so: $(OBJ)
	$(CC) $^ -shared -o $@ $(SOFLAGS)

python: pyIMSRG.so

pyIMSRG.so: $(OBJ) pyIMSRG.o
	$(CC) $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) -shared $(PYTHONFLAGS)

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) $(PYTHONFLAGS)

clean:
	rm -f *.o *.so
