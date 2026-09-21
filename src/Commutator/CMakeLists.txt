add_library(
    Commutator   OBJECT
    Commutator.hh Commutator.cc
    DaggerCommutators.hh DaggerCommutators.cc
    TensorCommutators.hh TensorCommutators.cc
    IMSRG3Commutators.hh IMSRG3Commutators.cc
    FactorizedDoubleCommutator.hh FactorizedDoubleCommutator.cc
)
target_include_directories(
    Commutator
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
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


