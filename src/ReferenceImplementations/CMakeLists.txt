
add_library(
    ReferenceImplementations  OBJECT
    ReferenceImplementations.hh 
    ReferenceImplementations.cc
    HFMBPT4.cc
    IMSRG2Commutators.cc
    IMSRG3Commutators.cc
    IMSRG2TensorCommutators.cc
    IMSRG3TensorCommutators.cc
    FactorizedNestedCommutators.cc
)

target_include_directories(
   ReferenceImplementations
   PUBLIC
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/..
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

