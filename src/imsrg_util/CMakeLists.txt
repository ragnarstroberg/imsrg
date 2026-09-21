
add_library(
    imsrg_util OBJECT
    imsrg_util.hh
    imsrg_util.cc
    GaussLaguerre.hh
    M0nu.hh
    M0nu.cc
    DarkMatterNREFT.hh
    DarkMatterNREFT.cc
    Pwd.hh
    Pwd.cc
    Atomic.cc
#    Helicity.hh Helicity.cc
)
target_include_directories(
    imsrg_util
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
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
