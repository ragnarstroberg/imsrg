CXX = g++
#CXX = clang++

# Compiler flags:
#  -DNO_ODE compiles without boost/ode package used for flow equation solver
#  -DNO_HDF5 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flag
.PHONY: version
ALL = version libIMSRG.so imsrg++
INSTDIR = $(HOME)
INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp -fPIC
LIBS = -lopenblas -lgsl -lz
#LIBS = -lblas -lgsl -lz

# I assume we're running on linux
OS = LINUX
# But in case we're crazy enough to run on MacOS, might as well check...
ifneq (,$(findstring arwin,$(shell uname)))
  OS = MACOS
endif


# Tag the build with the current version from git, so we know where all the awesome answers came from.
VERSIONCOMMAND = 
BUILDVERSION = $(shell git branch -v | grep '*' | awk '{printf "%s_%s",$$2,$$3}' )
ifneq (,$(findstring detached,$(BUILDVERSION)))
  BUILDVERSION = $(shell git branch -v | grep '*' | awk '{printf "HEAD_%s",$$6}'  )
endif
OLDBUILD = $(shell awk '{print $$3}' version.hh)

ifneq ($(OLDBUILD),"$(BUILDVERSION)")
  VERSIONCOMMAND = echo "\#define BUILDVERSION \"$(BUILDVERSION)\"" > version.hh
endif

#FLAGS += -DBUILDVERSION=\"$(BUILDVERSION)\"

ifeq ($(DEBUG),on)
 FLAGS     = -march=native -std=c++11 -fopenmp -fPIC -g
endif

#WITHHALF = 1

#ifeq (1,$(WITHHALF))
  INCLUDE += -I./half/include
#  FLAGS += -DWITHHALF=1
#endif

WARNFLAGS = -Wall -Wno-comment -Wno-deprecated-declarations -Wno-missing-braces
FLAGS    += $(WARNFLAGS)

SOFLAGS = $(FLAGS)

PYTHON_MAJOR_VERSION = $(shell python -c "import sys; sys.stdout.write('{}'.format(sys.version_info.major))")
ifneq (,$(shell which python3)) # If it's there, default to python3
  PYTHON_MAJOR_VERSION = 3
endif
PYTHON_CONFIG = python$(PYTHON_MAJOR_VERSION)-config
ifeq (,$(shell which $(PYTHON_CONFIG) 2> /dev/null))
  PYTHON_CONFIG = python-config
endif
PYTHON_CFLAGS = $(shell $(PYTHON_CONFIG) --cflags )
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags )
PYTHON_INCLUDE = -I./pybind11/include
PYTHONFLAGS =  $(PYTHON_CFLAGS) -Wl,--no-as-needed
PYTHON_COMMAND = python$(PYTHON_MAJOR_VERSION)

THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; elif [ `hostname|grep oak` ]; then echo oak; elif [ `hostname|grep cedar` ]; then echo cedar; elif [ `hostname|grep mox` ]; then echo mox; elif [ `hostname|grep bebop` ]; then echo bebop; else echo other; fi)



ifeq ($(OS),MACOS)
  FLAGS     = -Xpreprocessor -fopenmp -O3  -std=c++11 -fPIC
  LIBS += -lomp
  FLAGS += -DNO_x86  # don't load x86intrin.h header. I don't pretent to fully understand why this breaks things.
  PYTHONFLAGS =  $(shell $(PYTHON_CONFIG) --cflags | sed -e 's/-arch i386//')
  PYTHON_LDFLAGS += -undefined dynamic_lookup # This avoids trouble with missing symbols (see pybind11.readthedocs.io/en/stable/compiling.html )
  PYTHONFLAGS += -Wno-unused-value # suppress annoying warnings about stuff in the pybind11 library
  PYTHONFLAGS += -Wno-deprecated-declarations
endif


ifeq ($(HDF5),on)
 LIBS += -lhdf5_cpp
 NEWHDF5LOC = $(shell if [ -d /usr/include/hdf5/serial ]; then echo yes; else echo no; fi)
 ifeq ($(NEWHDF5LOC),yes)
   INCLUDE += -I/usr/include/hdf5/serial
   LIBS += -lhdf5_serial
 else
  LIBS += -lhdf5
 endif
else
 FLAGS += -DNO_HDF5   # By default, don't bother building with HDF5, since it's not used very often and can make building a pain
endif



ifeq ($(THEHOST),other)  # default options. assumes boost and python are set up nicely.
 LIBS += -llapack
# FLAGS += -DOPENBLAS_NOUSEOMP=1
# ifneq ($(PYTHON),off)
#  ALL += pyIMSRG.so
# endif
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
 PYTHONFLAGS := -L/opt/anaconda/lib $(PYTHONFLAGS)
endif

ifeq ($(THEHOST),cronos)
 FLAGS   += -DOLD_BOOST=1
 SOFLAGS += -DOLD_BOOST=1
 INCLUDE += -I$(HOME)/include
 LIBS    += -L$(HOME)/lib

endif

ifeq ($(THEHOST),oak)
  NEWLIBS := $(filter-out -lopenblas -lblas,$(LIBS))
  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
  LIBS +=  -lmkl_avx2  -lmkl_def
# ifneq ($(PYTHON),off)
#  ALL += pyIMSRG.so
#  endif
endif

ifeq ($(THEHOST),mox)
  NEWLIBS := $(filter-out -lopenblas -lblas,$(LIBS))
  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif



ifeq ($(THEHOST),cedar)
#  NEWLIBS := $(filter-out -lopenblas,$(LIBS))
#  NEWLIBS := $(filter-out -lblas,$(LIBS))
#  LIBS = $(NEWLIBS) -lopenblas
#  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
# ifneq ($(PYTHON),off)
#  ALL += pyIMSRG.so
#  endif
endif

ifeq ($(THEHOST),bebop)
  NEWLIBS := $(filter-out -lopenblas -lblas,$(LIBS))
  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
  LIBS +=  -lmkl_avx2  -lmkl_def
# ifneq ($(PYTHON),off)
#  ALL += pyIMSRG.so
#  endif
endif




ifneq ($(PYTHON),off)
 ALL += pyIMSRG.so
 endif


all: $(ALL)
	@echo Building with build version $(BUILDVERSION)

OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
      HartreeFock.o imsrg_util.o Generator.o IMSRGSolver.o AngMom.o AngMomCache.o\
      IMSRGProfiler.o Commutator.o HFMBPT.o  RPA.o\
      M0nu.o DarkMatterNREFT.o Jacobi3BME.o UnitTest.o \
      TwoBodyChannel.o ThreeBodyChannel.o \
      ThreeBodyStorage.o ThreeBodyStorage_pn.o ThreeBodyStorage_iso.o \
      ThreeBodyStorage_no2b.o ThreeBodyStorage_mono.o ThreeLegME.o  \
      boost_src/gzip.o boost_src/zlib.o
#      ThreeBodyStorage_no2b.o ThreeLegME.o ThreeBodyMENO2B.o \
#      ThreeBodyMEpn.o ThreeLegME.o ThreeBodyMENO2B.o \


boost_src/%.o: boost_src/%.cpp
	$(CXX) -c $^ -o $@ $(INCLUDE) $(FLAGS)


%.o: %.cc %.hh
	$(CXX) -c $*.cc -o $@ $(INCLUDE) $(FLAGS)

%.o: %.cc
	$(CXX) -c $*.cc -o $@ $(INCLUDE) $(FLAGS)


libIMSRG.so: $(OBJ)
	$(CXX) $^ -shared -o $@ $(SOFLAGS) $(LIBS)

python: pyIMSRG.so

pyIMSRG.so: $(OBJ)  pyIMSRG.o version.hh
	$(CXX) $(FLAGS) $(INCLUDE) $(PYTHON_INCUDE) $(PYTHONFLAGS) $(PYTHON_LDFLAGS) $(LIBS) -shared  $(OBJ) pyIMSRG.o -o $@

pyIMSRG.o: pyIMSRG.cc
	$(CXX) -c  $(FLAGS) $(INCLUDE) $(PYTHON_INCLUDE) $(PYTHONFLAGS)  $^ -o $@

imsrg++: imsrg++.cc libIMSRG.so Parameters.hh version.hh
	$(CXX) $(INCLUDE) $< -o $@ $(FLAGS) -L$(PWD) -lIMSRG $(LIBS)


version:
	$(VERSIONCOMMAND)
#	echo "#define BUILDVERSION \"$(BUILDVERSION)\"" > version.hh
#
clean:
	rm -f *.o *.so boost_src/*.o

test:
	$(PYTHON_COMMAND) -c "import pyIMSRG; ms=pyIMSRG.ModelSpace(2,'He4','He4'); ut=pyIMSRG.UnitTest(ms); passed=ut.SanityCheck()"

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
	ln -nsf $(PWD)/armadillo $(INSTDIR)/include/armadillo
	@printf "\n\nDone installing.\n\n"
	@echo '*********************************************************************'
	@echo '* Make sure libIMSRG.so is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
	@echo '*********************************************************************'



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
