CC = g++

# Compiler flags:
#  -DNO_ODE=1 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flat

INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp -fPIC -ffast-math
SOFLABS   = -O3 -march=native -std=c++11 -fopenmp -fPIC -ffast-math
WARNFLAGS = -Wall -Wno-comment
FLAGS    += $(WARNFLAGS)
SOFLAGS  += $(WARNFLAGS)

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

ifneq (`echo $(HOSTNAME) | grep jrl` , "") # specific options for jureca cluster
 LIBS = -lopenblas -lgsl -lz -lhdf5_cpp -lhdf5
 FLAGS += -DOPENBLAS_NOUSEOMP=1
 ALL = libIMSRG.so
else
ifneq (`echo $(HOSTNAME) | grep cougar` "") # specific options for cougar cluster
 LIBS = -lopenblas -llapack -lgsl -lz -lhdf5_cpp -lhdf5
 ALL = libIMSRG.so
else  # default options
 LIBS = -lopenblas -llapack -lgsl -lz -lhdf5_cpp -lhdf5
 ALL = pyIMSRG.so libIMSRG.so
endif
endif

all: $(ALL)
	
OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
      HartreeFock.o Generator.o IMSRGSolver.o AngMom.o imsrg_util.o\
      boost_src/gzip.o boost_src/zlib.o

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/gzip.o: boost_src/gzip.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/zlib.o: boost_src/zlib.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 


%.o: %.cc %.hh
	$(CC) -c $*.cc -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 



libIMSRG.so: $(OBJ)
	$(CC) $^ -shared -o $@ $(SOFLAGS)

python: pyIMSRG.so

pyIMSRG.so: $(OBJ)  pyIMSRG.o
	$(CC) $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) -shared $(PYTHONFLAGS)

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) $(PYTHONFLAGS)

clean:
	rm -f *.o *.so boost_src/*.o
