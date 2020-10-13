cmake_minimum_required(VERSION 3.7)
project(imsrg++)

execute_process( COMMAND git branch -v COMMAND grep \* COMMAND awk "{printf \"%s_%s\",$2,$3}"  OUTPUT_VARIABLE BUILDVERSION)
MESSAGE( "BUILDVERSION is ${BUILDVERSION}")

add_custom_target( version.hh ALL  COMMAND git branch -v | grep \* | awk "{printf \"#define BUILDVERSION \\\"%s_%s\\\"\\n\",$2,$3}" > version.hh VERBATIM USES_TERMINAL)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(INSTDIR $ENV{HOME})
set(CMAKE_INSTALL_MESSAGE LAZY)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -DBUILDVERSION=\\\"${BUILDVERSION}\\\"")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -DBUILDVERSION=\\\"${BUILDVERSION}\\\"")
#MESSAGE( "CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

set(splashscreen 
"                                                     ____                                     
             _________________          _____________/   /\\               _________________    
           /____/_____/_____/|         /____/_____/ /___/  \\             /____/_____/_____/|   
          /____/_____/__G_ /||        /____/_____/|/   /\\  /\\           /____/_____/____ /||   
         /____/_____/__+__/|||       /____/_____/|/ G /  \\/  \\         /____/_____/_____/|||   
        |     |     |     ||||      |     |     |/___/   /\\  /\\       |     |     |     ||||   
        |  I  |  M  |     ||/|      |  I  |  M  /   /\\  /  \\/  \\      |  I  |  M  |     ||/|   
        |_____|_____|_____|/||      |_____|____/ + /  \\/   /\\  /      |_____|_____|_____|/||   
        |     |     |     ||||      |     |   / __/   /\\  /  \\/       |     |     |     ||||   
        |  S  |  R  |     ||/|      |  S  |   \\   \\  /  \\/   /        |  S  |  R  |  G  ||/|   
        |_____|_____|_____|/||      |_____|____\\ __\\/   /\\  /         |_____|_____|_____|/||   
        |     |     |     ||||      |     |     \\   \\  /  \\/          |     |     |     ||||   
        |     |  +  |     ||/       |     |  +  |\\ __\\/   /           |     |  +  |  +  ||/    
        |_____|_____|_____|/        |_____|_____|/\\   \\  /            |_____|_____|_____|/     
                                                   \\___\\/                                      
                                                                                               
")

#message(${splashscreen})
add_executable(imsrg++ imsrg++.cc)


add_library(OBJ OBJECT  ModelSpace.cc TwoBodyME.cc ThreeBodyME.cc Operator.cc  ReadWrite.cc
                        HartreeFock.cc imsrg_util.cc Generator.cc IMSRGSolver.cc AngMom.cc
                        IMSRGProfiler.cc Commutator.cc HFMBPT.cc
                        M0nu.cc DarkMatterNREFT.cc Jacobi3BME.cc UnitTest.cc 
                        TwoBodyChannel.cc ThreeBodyChannel.cc 
                        ThreeBodyStorage.cc ThreeBodyStorage_pn.cc ThreeBodyStorage_iso.cc 
                        ThreeBodyStorage_no2b.cc ThreeLegME.cc
                        boost_src/gzip.cpp boost_src/zlib.cpp)
add_library(IMSRG SHARED)


target_link_libraries(imsrg++ PUBLIC IMSRG)
target_include_directories(OBJ PUBLIC armadillo)
target_include_directories(OBJ PUBLIC half/include)
target_link_libraries(IMSRG PUBLIC OBJ)
target_link_libraries(imsrg++ PUBLIC IMSRG)
#Links to the different needed packages
find_package(GSL REQUIRED)
if (GSL_FOUND) 
  target_link_libraries(OBJ PUBLIC GSL::gsl GSL::gslcblas)
endif(GSL_FOUND)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OBJ PUBLIC OpenMP::OpenMP_CXX)
endif(OpenMP_CXX_FOUND)

find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIR})
#  set(_hdf5_libs hdf5 hdf5_cpp)
  set(_hdf5_libs ${HDF5_LIBRARIES})
  target_link_libraries(OBJ PUBLIC ${_hdf5_libs})
else()
  target_compile_options(OBJ PUBLIC -DNO_HDF5)
endif(HDF5_FOUND)

#set(BLA_VENDOR OpenBLAS)
#set(BLA_PREFER_PKGCONFIG)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(BLAS_FOUND AND LAPACK_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  target_link_libraries(OBJ PUBLIC ${lapackblas_libraries})
endif(BLAS_FOUND AND LAPACK_FOUND)
MESSAGE("labpack/blas libraries:  ${lapackblas_libraries}")

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  target_link_libraries(OBJ PUBLIC ZLIB::ZLIB)
endif(ZLIB_FOUND)
find_package(PythonInterp)
find_package(PythonLibs)


####Installation command###
install(CODE "message(\"${splashscreen}\")")
# Install imsrg executable
install(TARGETS imsrg++ 
        DESTINATION ${INSTDIR}/bin 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE) 

# Install library
install(TARGETS IMSRG  
        DESTINATION ${INSTDIR}/lib 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Install the header files
file(GLOB HEADERS *.hh)
install(FILES ${HEADERS} 
        DESTINATION ${INSTDIR}/include
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(DIRECTORY armadillo 
        DESTINATION ${INSTDIR}/include 
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(CODE "MESSAGE(\"Done installing.\")")
install(CODE "MESSAGE(\"*********************************************************************\")")
install(CODE "MESSAGE(\"* Make sure libIMSRG.so is in your LIBRARY_PATH and LD_LIBRARY_PATH *\")")
install(CODE "MESSAGE(\"*********************************************************************\")")



###Compile pyIMSRG shared library if it finds python
if (PythonLibs_FOUND AND PythonInterp_FOUND)
  add_subdirectory(pybind11)
  pybind11_add_module(pyIMSRG pyIMSRG.cc)
  target_link_libraries(pyIMSRG PUBLIC OBJ)
  target_include_directories(pyIMSRG PUBLIC armadillo)
  target_include_directories(pyIMSRG PUBLIC half/include)
  set_target_properties(pyIMSRG PROPERTIES CXX_VISIBILITY_PRESET default)
  install(TARGETS pyIMSRG  DESTINATION ${INSTDIR}/lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ)
endif(PythonLibs_FOUND AND PythonInterp_FOUND)


add_custom_target(test COMMAND python -c "import pyIMSRG; ms=pyIMSRG.ModelSpace(2,'He4','He4'); ut=pyIMSRG.UnitTest(ms); passed=ut.SanityCheck()" VERBATIM USES_TERMINAL)
