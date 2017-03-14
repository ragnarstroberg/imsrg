CC = g++

# Compiler flags:
#  -DNO_ODE=1 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flag

INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp -fPIC 
SOFLAGS   = -O3 -march=native -std=c++11 -fopenmp -fPIC 
INSTDIR = $(HOME)

BUILDVERSION = $(shell git branch -v | grep '*' | awk '{printf "%s_%s",$$2,$$3}' )
FLAGS += -DBUILDVERSION=\"$(BUILDVERSION)\"

ifeq ($(DEBUG),on)
 FLAGS     = -march=native -std=c++11 -fopenmp -fPIC -g 
 SOFLAGS   = -march=native -std=c++11 -fopenmp -fPIC -g
endif

WARNFLAGS = -Wall -Wno-comment -Wno-deprecated-declarations
FLAGS    += $(WARNFLAGS)
SOFLAGS  += $(WARNFLAGS)

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; else echo other; fi)

ALL = libIMSRG.so imsrg++
LIBS = -lopenblas -lgsl -lz -lhdf5_cpp -lhdf5

ifeq ($(THEHOST),jureca) # specific options for jureca cluster
# FLAGS += -DOPENBLAS_NOUSEOMP=1
 SOFLAGS += -fuse-ld=bfd
endif
ifeq ($(THEHOST),cougar) # specific options for cougar cluster
 LIBS += -llapack
 PYTHONFLAGS = -I/opt/anaconda/include/python2.7 -L/opt/boost/1.58.0/lib -lboost_python #-lpython 
endif
ifeq ($(THEHOST),cronos)
 FLAGS += -DOLD_BOOST=1
 SOFLAGS += -DOLD_BOOST=1
 INCLUDE += -I$(HOME)/include
 LIBS += -L$(HOME)/lib
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
  endif
endif
ifeq ($(THEHOST),other)  # default options. assumes boost and python are set up nicely.
 LIBS += -llapack
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
 endif
endif

# The new version of Ubuntu puts the hdf5 libraries in a different location
NEWHDF5LOC = $(shell if [ -d /usr/include/hdf5/serial ]; then echo yes; else echo no; fi)
ifeq ($(NEWHDF5LOC),yes)
  INCLUDE += -I/usr/include/hdf5/serial
  NEWLIBS := $(filter-out -lhdf5,$(LIBS))
  LIBS = $(NEWLIBS) -lhdf5_serial
endif

all: $(ALL)
	@echo Building with build version $(BUILDVERSION)

install: splashscreen
	@if [ ! -d $(INSTDIR)/lib ] ; then \
	  mkdir $(INSTDIR)/lib; \
	fi
	@if [ ! -d $(INSTDIR)/include ] ; then \
	  mkdir $(INSTDIR)/include; \
	fi 
	ln -sf $(PWD)/libIMSRG.so $(INSTDIR)/lib/libIMSRG.so
	ln -sfT $(PWD)/armadillo $(INSTDIR)/include/armadillo
	@if [ -f pyIMSRG.so ] ; then\
	  ln -sf $(PWD)/pyIMSRG.so $(INSTDIR)/lib/pyIMSRG.so;\
	fi
	@for x in *.hh; do \
	 echo linking $(PWD)/$$x  '=>'  $(INSTDIR)/include/$$x;\
	 ln -sf $(PWD)/$$x $(INSTDIR)/include/$$x; \
	done
	@if [ -d $(INSTDIR)/bin ] ; then \
	  echo linking $(PWD)/imsrg++  '=>' $(INSTRDIR)/bin/imsrg++;\
	  ln -sf $(PWD)/imsrg++ $(INSTDIR)/bin/imsrg++;\
	fi
	@printf "\n\nDone installing.\n\n"
	@echo '*********************************************************************'
	@echo '* Make sure libIMSRG.so is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
	@echo '*********************************************************************'
	 

OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
      HartreeFock.o imsrg_util.o Generator.o IMSRGSolver.o AngMom.o\
      boost_src/gzip.o boost_src/zlib.o  IMSRGProfiler.o 

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

imsrg++: imsrg++.cc libIMSRG.so Parameters.hh
	$(CC) $(INCLUDE) $< -o $@ $(FLAGS) -L$(PWD) -lIMSRG $(LIBS)

clean:
	rm -f *.o *.so boost_src/*.o



splashscreen:
	@printf "                                                      ____                                     \n"   
	@printf "             _________________          _____________/   /\               _________________    \n"   
	@printf "           /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|   \n"  
	@printf "          /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||   \n"  
	@printf "         /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||   \n"  
	@printf "        |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||   \n"  
	@printf "        |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|   \n"  
	@printf "        |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||   \n"  
	@printf "        |     |     |     ||||      |     |   / __/   /\  /  \/       |     |     |     ||||   \n"  
	@printf "        |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|   \n"  
	@printf "        |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||   \n"  
	@printf "        |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||   \n"  
	@printf "        |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/    \n"  
	@printf "        |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/     \n"  
	@printf "                                                   \___\/                                      \n"  
	@printf "                                                                                               \n"  
