CC = g++
#CC = clang++

# Compiler flags:
#  -DNO_ODE compiles without boost/ode package used for flow equation solver
#  -DNO_HDF5 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flag

ALL = libIMSRG.so imsrg++
INSTDIR = $(HOME)
INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp -fPIC 
LIBS = -lopenblas -lgsl -lz

# I assume we're running on linux
OS = LINUX
# But in case we're crazy enough to run on MacOS, might as well check...
ifneq (,$(findstring arwin,$(shell uname))) 
  OS = MACOS
endif


# Tag the build with the current version from git, so we know where all the awesome answers came from.
BUILDVERSION = $(shell git branch -v | grep '*' | awk '{printf "%s_%s",$$2,$$3}' )
ifneq (,$(findstring detached,$(BUILDVERSION)))
  BUILDVERSION = $(shell git branch -v | grep '*' | awk '{printf "HEAD_%s",$$6}'  )
endif

ifeq ($(DEBUG),on)
 FLAGS     = -march=native -std=c++11 -fopenmp -fPIC -g 
endif

FLAGS += -DBUILDVERSION=\"$(BUILDVERSION)\"

WARNFLAGS = -Wall -Wno-comment -Wno-deprecated-declarations -Wno-missing-braces
FLAGS    += $(WARNFLAGS)

SOFLAGS = $(FLAGS)

#PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 
PYTHONFLAGS = -I./pybind11/include $(shell python-config --cflags ) -Wl,--no-as-needed
PYTHONLIBS = -I./pybind11/include $(shell python-config  --ldflags) 

#THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; elif [ `hostname|grep oak` ]; then echo oak; else echo other; fi)
THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; elif [ `hostname|grep oak` ]; then echo oak; elif [ `hostname|grep cedar` ]; then echo cedar; else echo other; fi)



ifeq ($(OS),MACOS)
#  FLAGS     = -Xpreprocessor -fopenmp -O3 -march=native -std=c++11 -fPIC 
  FLAGS     = -Xpreprocessor -fopenmp -O3  -std=c++11 -fPIC 
  LIBS += -lomp
  PYTHONFLAGS = -I./pybind11/include $(shell python-config --cflags | sed -e 's/-arch i386//')
endif


ifeq ($(HDF5),on)
 LIBS += -lhdf5_cpp -lhdf5
else
 FLAGS += -DNO_HDF5   # By default, don't bother building with HDF5, since it's not used very often and can make building a pain
endif



ifeq ($(THEHOST),jureca) # specific options for jureca cluster
# FLAGS += -DOPENBLAS_NOUSEOMP=1
 SOFLAGS += -fuse-ld=bfd
 PYTHONFLAGS := $(filter-out -ftz,$(PYTHONFLAGS))
 PYTHONFLAGS := $(filter-out -fp-speculation=safe,$(PYTHONFLAGS))
 PYTHONFLAGS := $(filter-out -fp-model,$(PYTHONFLAGS))
 PYTHONFLAGS := $(filter-out source,$(PYTHONFLAGS))
 PYTHONFLAGS := $(filter-out -xHost,$(PYTHONFLAGS))
endif
ifeq ($(THEHOST),cougar) # specific options for cougar cluster
 LIBS += -llapack
# PYTHONFLAGS = -I/opt/anaconda/include/python2.7 -L/opt/boost/1.58.0/lib -lboost_python #-lpython 
  PYTHONFLAGS := -L/opt/anaconda/lib $(PYTHONFLAGS)
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

ifeq ($(THEHOST),oak)
#  FLAGS += -DOLD_BOOST
#  FLAGS += -DOPENBLAS_NOUSEOMP=1
  NEWLIBS := $(filter-out -lopenblas,$(LIBS))
  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
  endif
endif

ifeq ($(THEHOST),cedar)
#  NEWLIBS := $(filter-out -lopenblas,$(LIBS))
#  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
 LIBS += -lmkl -lmkl_gnu_thread -lmkl_core -lpthread -lgomp -liomp5
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
  endif
endif



ifeq ($(HDF5),on)
# The new version of Ubuntu puts the hdf5 libraries in a different location
 NEWHDF5LOC = $(shell if [ -d /usr/include/hdf5/serial ]; then echo yes; else echo no; fi)
 ifeq ($(NEWHDF5LOC),yes)
   INCLUDE += -I/usr/include/hdf5/serial
   NEWLIBS := $(filter-out -lhdf5,$(LIBS))
   LIBS = $(NEWLIBS) -lhdf5_serial
 endif
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
	@if [ ! -d $(INSTDIR)/bin ] ; then \
	  mkdir $(INSTDIR)/bin; \
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
      IMSRGProfiler.o \
      boost_src/gzip.o boost_src/zlib.o  

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/gzip.o: boost_src/gzip.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(FLAGS) 
#	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/zlib.o: boost_src/zlib.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(FLAGS) 
#	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 


%.o: %.cc %.hh
	$(CC) -c $*.cc -o $@ $(INCLUDE) $(FLAGS) 
#	$(CC) -c $*.cc -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 



libIMSRG.so: $(OBJ)
	$(CC) $^ -shared -o $@ $(SOFLAGS) $(LIBS)

python: pyIMSRG.so

pyIMSRG.so: $(OBJ)  pyIMSRG.o
	$(CC) $(FLAGS) $(INCLUDE) $(PYTHONFLAGS) $(PYTHONLIBS) $(LIBS) -shared  $^ -o $@

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c  $(FLAGS) $(INCLUDE) $(PYTHONFLAGS)  $^ -o $@
#	$(CC) -c  $(FLAGS) $(INCLUDE) $(PYTHONFLAGS) $(LIBS)  $^ -o $@

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
