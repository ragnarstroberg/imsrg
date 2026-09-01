cmake_minimum_required(VERSION 3.18)
set(INSTDIR $ENV{HOME})
set(BLA_VENDOR OpenBLAS)
# set(CMAKE_SKIP_INSTALL_RPATH ON)
# set(CMAKE_SYSTEM_PREFIX_PATH $ENV{EBROOTGENTOO})
set(CMAKE_CURRENT_BINARY_DIR "../build")

option(IMSRG_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

# Place to put any settings for quirks of your specific environment
# e.g. maybe you want a different BLAS, or maybe you need to
# explicitly set some additional variables for whatever reason
include(local_settings.cmake OPTIONAL)

# Determine if imsrg++ is built as a subproject (using add_subdirectory)
# or if it is the main project.
if(NOT DEFINED IMSRGPLUSPLUS_MAIN_PROJECT)
  set(IMSRGPLUSPLUS_MAIN_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(IMSRGPLUSPLUS_MAIN_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

# Set the default CMAKE_BUILD_TYPE to Release.
# This should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if(IMSRGPLUSPLUS_MAIN_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(imsrg++ VERSION 0.1.0)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(ZLIB REQUIRED)
# find_package(Armadillo REQUIRED)
# include_directories(${ARMADILLO_INCLUDE_DIRS})


add_compile_options(-frounding-math)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(IMSRGPLUSPLUS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(IMSRG_BUILD_PYTHON_BINDINGS)
    # We only need this for the Python bindings.
    # In principle, everything should be faster if we could avoid this.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif()

add_subdirectory(extern/armadillo)
add_subdirectory(extern/half)

add_subdirectory(extern/boost_src)

# add_subdirectory(profiling)


if(IMSRG_BUILD_PYTHON_BINDINGS)
    add_subdirectory(extern/pybind11)
endif()

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cc
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cc
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/git_version.cmake)


add_library(
    AngMom  OBJECT
    AngMom.hh AngMom.cc
    AngMomCache.hh AngMomCache.cc
)
target_link_libraries(
    AngMom
    PRIVATE
    GSL::gsl
)


add_library(
    IMSRGProfiler OBJECT
    IMSRGProfiler.hh IMSRGProfiler.cc
)
target_link_libraries(
    IMSRGProfiler
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    ModelSpace  OBJECT
    ModelSpace.hh ModelSpace.cc
    TwoBodyChannel.cc
    ThreeBodyChannel.cc
)
target_link_libraries(
    ModelSpace
    PUBLIC
    armadillo
    IMSRGProfiler
    AngMom
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    PhysicalConstants
    INTERFACE
    PhysicalConstants.hh
)

add_library(
    TwoBodyME  OBJECT
    TwoBodyME.hh TwoBodyME.cc
)
target_link_libraries(
    TwoBodyME
    PUBLIC
    ModelSpace
    PRIVATE
    AngMom
    PhysicalConstants
)

add_library(
    ThreeLegME OBJECT
    ThreeLegME.hh ThreeLegME.cc
)
target_link_libraries(
    ThreeLegME
    PUBLIC
    armadillo
    ModelSpace
    PRIVATE
    # ${ARMADILLO_LIBRARIES}
    PhysicalConstants
)


add_library(
    ThreeBodyME OBJECT
    ThreeBodyME.hh ThreeBodyME.cc
    ThreeBodyStorage.hh ThreeBodyStorage.cc
    ThreeBodyStorage_pn.hh ThreeBodyStorage_pn.cc
    ThreeBodyStorage_no2b.hh ThreeBodyStorage_no2b.cc
    ThreeBodyStorage_mono.hh ThreeBodyStorage_mono.cc
    ThreeBodyStorage_iso.hh ThreeBodyStorage_iso.cc
)
target_link_libraries(
    ThreeBodyME
    PUBLIC
    ModelSpace
    IMSRGBoostZip
    HalfPrecision
    PRIVATE
    IMSRGProfiler
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    Parameters
    INTERFACE
)

add_library(
    Operator  OBJECT
    Operator.hh Operator.cc
)
target_link_libraries(
    Operator
    PUBLIC
    armadillo
    ModelSpace
    TwoBodyME
    ThreeBodyME
    ThreeLegME
    IMSRGProfiler
    PRIVATE
    GSL::gsl
    OpenMP::OpenMP_CXX
    PhysicalConstants
    AngMom
)

add_library(
    Version  OBJECT
    version.hh
    ${CMAKE_CURRENT_BINARY_DIR}/version.cc
)
target_include_directories(
    Version
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(
    ReadWrite  OBJECT
    ReadWrite.hh ReadWrite.cc
)
target_compile_definitions(
    ReadWrite
    PUBLIC
    NO_HDF5="true"
)
target_link_libraries(
    ReadWrite
    PUBLIC
    ModelSpace
    Operator
    Jacobi3BME
    ThreeBodyME
    PRIVATE
    OpenMP::OpenMP_CXX
    IMSRGBoostZip
    AngMom
    PhysicalConstants
    Version
)

add_library(
    BCH  OBJECT
    BCH.hh BCH.cc
)
target_link_libraries(
    BCH
    PUBLIC
    Operator
    PRIVATE
    Commutator
    OpenMP::OpenMP_CXX
)

add_library(
    IMSRGSolver  OBJECT
    IMSRGSolver.hh IMSRGSolver.cc
    IMSRGSolverPV.hh IMSRGSolverPV.cc
)
target_link_libraries(
    IMSRGSolver
    PUBLIC
    Operator
    Generator
    IMSRGProfiler
    ReadWrite
    PRIVATE
    BCH
    Commutator
    Boost::boost
)

#add_library(
#    IMSRGSolverPV  OBJECT
#    IMSRGSolverPV.hh IMSRGSolverPV.cc
#)
#target_link_libraries(
#    IMSRGSolverPV
#    PUBLIC
#    IMSRGSolver
#    GeneratorPV
#    PRIVATE
#    Commutator
#)

add_library(
    Commutator   OBJECT
    Commutator.hh Commutator.cc
    DaggerCommutators.hh DaggerCommutators.cc
    TensorCommutators.hh TensorCommutators.cc
    IMSRG3Commutators.hh IMSRG3Commutators.cc
    FactorizedDoubleCommutator.hh FactorizedDoubleCommutator.cc
)
target_link_libraries(
    Commutator
    PUBLIC
    armadillo
    Operator
    TwoBodyME
    ThreeLegME
    PRIVATE
    OpenMP::OpenMP_CXX
    # ${ARMADILLO_LIBRARIES}
    PhysicalConstants
    AngMom
)


add_library(
    Jacobi3BME OBJECT
    Jacobi3BME.hh Jacobi3BME.cc
)
target_link_libraries(
    Jacobi3BME
    PUBLIC
    ModelSpace
    PRIVATE
    AngMom
    IMSRGProfiler
    OpenMP::OpenMP_CXX
)

add_library(
    HartreeFock OBJECT
    HartreeFock.hh HartreeFock.cc
)
target_link_libraries(
    HartreeFock
    PUBLIC
    armadillo
    ModelSpace
    Operator
    IMSRGProfiler
    Jacobi3BME
    ThreeBodyME
    PRIVATE
    OpenMP::OpenMP_CXX
    GSL::gsl
    AngMom
    PhysicalConstants
)




add_library(
    imsrg_util OBJECT
    imsrg_util.hh imsrg_util.cc
    GaussLaguerre.hh
    M0nu.hh M0nu.cc
    DarkMatterNREFT.hh DarkMatterNREFT.cc
    Pwd.hh Pwd.cc
#    Helicity.hh Helicity.cc
)
target_link_libraries(
    imsrg_util
    PUBLIC
    ModelSpace
    Operator
    HartreeFock
    PhysicalConstants
    GSL::gsl
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    Generator OBJECT
    Generator.hh Generator.cc
    GeneratorPV.hh GeneratorPV.cc
)
target_link_libraries(
    Generator
    PUBLIC
    ModelSpace
    Operator
    PRIVATE
    GSL::gsl
    imsrg_util
    AngMom
#    Commutator
    PhysicalConstants
    OpenMP::OpenMP_CXX
)

#add_library(
#    GeneratorPV OBJECT
#    GeneratorPV.hh GeneratorPV.cc
#)
#target_link_libraries(
#    GeneratorPV
#    PUBLIC
#    ModelSpace
#    Operator
#    PRIVATE
#    Generator
#    GSL::gsl
#    imsrg_util
#    AngMom
##    Commutator
#    PhysicalConstants
#    OpenMP::OpenMP_CXX
#)

add_library(
    HFMBPT OBJECT
    HFMBPT.hh HFMBPT.cc
)
target_link_libraries(
    HFMBPT
    PUBLIC
    armadillo
    HartreeFock
    PRIVATE
    PhysicalConstants
    OpenMP::OpenMP_CXX
    Commutator
)



add_library(
    RPA OBJECT
    RPA.hh RPA.cc
)
target_link_libraries(
    RPA
    PUBLIC
    armadillo
    ModelSpace
    Operator
    PRIVATE
    PhysicalConstants
    AngMom
)

add_library(
    ReferenceImplementations  OBJECT
    ReferenceImplementations.hh ReferenceImplementations.cc
)
target_link_libraries(
    ReferenceImplementations
    PUBLIC
    ModelSpace
    Operator
#    Commutator
    PRIVATE
    armadillo
    PhysicalConstants
    Version
    imsrg_util
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    UnitTest  OBJECT
    UnitTest.hh UnitTest.cc
)
target_link_libraries(
    UnitTest
    PUBLIC
    ModelSpace
    Operator
    PRIVATE
    Commutator
    armadillo
    PhysicalConstants
    Version
    imsrg_util
    AngMom
    OpenMP::OpenMP_CXX
    ReferenceImplementations
)

add_library(
    IMSRG
    SHARED
)
target_link_libraries(
    IMSRG
#    INTERFACE
    Version
    ReadWrite
    ModelSpace
    TwoBodyME
    ThreeBodyME
    ThreeLegME
    Operator
    HartreeFock
    Commutator
    BCH
    IMSRGSolver
#    IMSRGSolverPV
    Generator
#    GeneratorPV
    imsrg_util
#    M0nu
#    Pwd
    AngMom
    IMSRGProfiler
    Jacobi3BME
#    DarkMatterNREFT
    HFMBPT
    UnitTest
    ReferenceImplementations
    PhysicalConstants
    RPA
)
set_target_properties(IMSRG PROPERTIES SUFFIX ".so")
set_target_properties(IMSRG PROPERTIES LINKER_LANGUAGE CXX)


add_executable(
    imsrg++
    imsrg++.cc
)
target_link_libraries(
    imsrg++
    PUBLIC
    IMSRG
    Parameters
    PhysicalConstants
    Version
    OpenMP::OpenMP_CXX
)


#target_link_libraries(
#    UnitTest
#    PUBLIC
#    IMSRG
#    Parameters
#    PhysicalConstants
#    Version
#    OpenMP::OpenMP_CXX
#)

if(IMSRG_BUILD_PYTHON_BINDINGS)
    pybind11_add_module(pyIMSRG pyIMSRG.cc)
    target_link_libraries(
        pyIMSRG
        PRIVATE
        IMSRG
    )
    set_target_properties(pyIMSRG PROPERTIES SUFFIX ".so")
    set_target_properties(pyIMSRG PROPERTIES POSITION_INDEPENDENT_CODE ON)
    add_custom_command(TARGET pyIMSRG PRE_LINK
		       COMMAND ${CMAKE_COMMAND} -E remove -f $<TARGET_FILE:pyIMSRG>
                       COMMENT "Unlinking old shared library to prevent SIGBUS on running jobs" )
endif()

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

include(CTest)
list(APPEND CMAKE_CTEST_ARGUMENTS "--verbose")
# Unit Tests
add_test(NAME Commutators2 COMMAND  python3 ../test/UnitTestCommutators2.py )
add_test(NAME Commutators3 COMMAND  python3 ../test/UnitTestCommutators3.py )
add_test(NAME Factorized   COMMAND  python3 ../test/UnitTestFactorized.py )
add_test(NAME Triples      COMMAND  python3 ../test/UnitTestTriples.py )
add_test(NAME Tensor       COMMAND  python3 ../test/UnitTestTensor.py )
add_test(NAME NormalOrdering    COMMAND  python3 ../test/UnitTestNormalOrdering.py )
