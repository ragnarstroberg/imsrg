########################################################################################################################
### IMSRG LIBRARY MAIN
########################################################################################################################

cmake_minimum_required(VERSION 3.22)

########################################################################################################################
### PROJECT DESCRIPTION AND META INFO
########################################################################################################################

set(IMSRG_PROJECT_NAME        "imsrg++")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(IMSRG_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(USE_ARMADILLO_EXTERN "Use version of armadillo in extern" ON)
########################################################################################################################
### DECLARE PROJECT
########################################################################################################################

project(imsrg++ VERSION 0.1.0)

set(IMSRG_SRC_DIR ${PROJECT_SOURCE_DIR}/src)


find_package(OpenMP REQUIRED)
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED)
if (NOT USE_ARMADILLO_EXTERN)
    find_package(Armadillo)
endif()

if(IMSRG_BUILD_PYTHON_BINDINGS)
    # We only need this for the Python bindings.
    # In principle, everything should be faster if we could avoid this.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif()

add_subdirectory(src)
add_subdirectory(extern)
if(IMSRG_BUILD_PYTHON_BINDINGS)
    pybind11_add_module(pyIMSRG src/pyIMSRG.cc)
    target_link_libraries(
        pyIMSRG
        PRIVATE
        IMSRG
    )
    #ADD tests to be run directly in the build process
    enable_testing()
    add_test (
        NAME UnitTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/UnitTest.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestCommutators2
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestCommutators2.py   
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test(
        NAME UnitTestCommutators3
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestCommutators3.py
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}   
        
    )
    add_test (
        NAME UnitTestFactorized
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestFactorized.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestTriples
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestTriples.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestTensor
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestTensor.py  
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()

