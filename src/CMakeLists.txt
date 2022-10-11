cmake_minimum_required(VERSION 3.18)

option(IMSRG_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

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


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(IMSRGPLUSPLUS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(IMSRG_BUILD_PYTHON_BINDINGS)
    # We only need this for the Python bindings.
    # In principle, everything should be faster if we could avoid this.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif()

add_subdirectory(armadillo)
add_subdirectory(half)

add_subdirectory(boost_src)

add_subdirectory(profiling)


if(IMSRG_BUILD_PYTHON_BINDINGS)
    add_subdirectory(pybind11)
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
    AngMom
    AngMom.hh AngMom.cc
)
target_link_libraries(
    AngMom
    PRIVATE
    GSL::gsl
)

add_library(
    AngMomCache
    AngMomCache.hh AngMomCache.cc
)
target_link_libraries(
    AngMomCache
    PRIVATE
    AngMom
)

add_library(
    IMSRGProfiler
    IMSRGProfiler.hh IMSRGProfiler.cc
)
target_link_libraries(
    IMSRGProfiler
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    ModelSpace
    ModelSpace.hh ModelSpace.cc
    TwoBodyChannel.cc
    ThreeBodyChannel.cc
)
target_link_libraries(
    ModelSpace
    PUBLIC
    armadillo
    IMSRGProfiler
    AngMomCache
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    PhysicalConstants
    INTERFACE
)

add_library(
    TwoBodyME
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
    ThreeLegME
    ThreeLegME.hh ThreeLegME.cc
)
target_link_libraries(
    ThreeLegME
    PUBLIC
    armadillo
    ModelSpace
    PRIVATE
    PhysicalConstants
)

add_library(
    ThreeBodyStorage
    ThreeBodyStorage.hh ThreeBodyStorage.cc
)
target_link_libraries(
    ThreeBodyStorage
    PUBLIC
    ModelSpace
)

add_library(
    ThreeBodyStorage_pn
    ThreeBodyStorage_pn.hh ThreeBodyStorage_pn.cc
)
target_link_libraries(
    ThreeBodyStorage_pn
    PUBLIC
    ThreeBodyStorage
    PRIVATE
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    ThreeBodyStorage_no2b
    ThreeBodyStorage_no2b.hh ThreeBodyStorage_no2b.cc
)
target_link_libraries(
    ThreeBodyStorage_no2b
    PUBLIC
    IMSRGBoostZip
    HalfPrecision
    ModelSpace
    ThreeBodyStorage
    PRIVATE
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    ThreeBodyStorage_mono
    ThreeBodyStorage_mono.hh ThreeBodyStorage_mono.cc
)
target_link_libraries(
    ThreeBodyStorage_mono
    PUBLIC
    IMSRGBoostZip
    HalfPrecision
    ModelSpace
    ThreeBodyStorage
    PRIVATE
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    ThreeBodyStorage_iso
    ThreeBodyStorage_iso.hh ThreeBodyStorage_iso.cc
)
target_link_libraries(
    ThreeBodyStorage_iso
    PUBLIC
    ThreeBodyStorage
    PRIVATE
    AngMom
)

add_library(
    ThreeBodyME
    ThreeBodyME.hh ThreeBodyME.cc
)
target_link_libraries(
    ThreeBodyME
    PUBLIC
    ModelSpace
    ThreeBodyStorage
    ThreeBodyStorage_iso
    ThreeBodyStorage_pn
    ThreeBodyStorage_no2b
    ThreeBodyStorage_mono
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
    Operator
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
    Version
    version.hh
    ${CMAKE_CURRENT_BINARY_DIR}/version.cc
)
target_include_directories(
    Version
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(
    ReadWrite
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
    IMSRGSolver
    IMSRGSolver.hh IMSRGSolver.cc
)
target_link_libraries(
    IMSRGSolver
    PUBLIC
    Operator
    Generator
    IMSRGProfiler
    ReadWrite
    PRIVATE
    Commutator
    Boost::boost
)

add_library(
    Commutator232
    Commutator232.hh Commutator232.cc
)
target_link_libraries(
    Commutator232
    PUBLIC
    armadillo
    Operator
    TwoBodyME
    ThreeLegME
    PRIVATE
    OpenMP::OpenMP_CXX
    PhysicalConstants
    AngMom
)

add_library(
    Commutator
    Commutator.hh Commutator.cc
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
    PhysicalConstants
    Commutator232
    AngMom
)

add_library(
    DarkMatterNREFT
    DarkMatterNREFT.hh DarkMatterNREFT.cc
)
target_link_libraries(
    DarkMatterNREFT
    PUBLIC
    ModelSpace
    Operator
    PRIVATE
    GSL::gsl
    PhysicalConstants
)

add_library(
    Jacobi3BME
    Jacobi3BME.hh Jacobi3BME.cc
)
target_link_libraries(
    Jacobi3BME
    PUBLIC
    ModelSpace
    PRIVATE
    AngMom
    # This is safe because HartreeFock is a private dep of Jacobi3BME.
    # If both deps were public, we would have a problem.
    HartreeFock
    IMSRGProfiler
    OpenMP::OpenMP_CXX
)

add_library(
    HartreeFock
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
    GaussLaguerre
    INTERFACE
)

add_library(
    M0nu
    M0nu.hh M0nu.cc
)
target_link_libraries(
    M0nu
    PUBLIC
    OpenMP::OpenMP_CXX
    GSL::gsl
    ModelSpace
    Operator
    imsrg_util
    PhysicalConstants
    PRIVATE
    AngMom
)

add_library(
    imsrg_util
    imsrg_util.hh imsrg_util.cc
)
target_link_libraries(
    imsrg_util
    PUBLIC
    ModelSpace
    Operator
    HartreeFock
    IMSRGSolver
    PhysicalConstants
    GSL::gsl
    PRIVATE
    OpenMP::OpenMP_CXX
    Commutator
    DarkMatterNREFT
    GaussLaguerre
    M0nu
)

add_library(
    Generator
    Generator.hh Generator.cc
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
    Commutator
    PhysicalConstants
    OpenMP::OpenMP_CXX
)

add_library(
    HFMBPT
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
)

add_library(
    RPA
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
    ReferenceImplementations
    ReferenceImplementations.hh ReferenceImplementations.cc
)
target_link_libraries(
    ReferenceImplementations
    PUBLIC
    ModelSpace
    Operator
    Commutator
    PRIVATE
    armadillo
    PhysicalConstants
    Version
    imsrg_util
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    UnitTest
    UnitTest.hh UnitTest.cc
)
target_link_libraries(
    UnitTest
    PUBLIC
    ModelSpace
    Operator
    Commutator
    PRIVATE
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
    INTERFACE
)
target_link_libraries(
    IMSRG
    INTERFACE
    ReadWrite
    ModelSpace
    TwoBodyME
    ThreeBodyME
    Operator
    HartreeFock
    Commutator
    IMSRGSolver
    Generator
    imsrg_util
    AngMom
    IMSRGProfiler
    Jacobi3BME
    DarkMatterNREFT
    HFMBPT
    UnitTest
    PhysicalConstants
    RPA
)

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

add_executable(
    unit_test
    unit_test.cc
)
target_link_libraries(
    unit_test
    PUBLIC
    IMSRG
    Parameters
    PhysicalConstants
    Version
    OpenMP::OpenMP_CXX
)

if(IMSRG_BUILD_PYTHON_BINDINGS)
    pybind11_add_module(pyIMSRG pyIMSRG.cc)
    target_link_libraries(
        pyIMSRG
        PRIVATE
        IMSRG
    )
endif()
