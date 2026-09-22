add_library(
    ThreeBodyME OBJECT
    ThreeBodyME.hh           ThreeBodyME.cc
    ThreeBodyStorage.hh      ThreeBodyStorage.cc
    ThreeBodyStorage_pn.hh   ThreeBodyStorage_pn.cc
    ThreeBodyStorage_no2b.hh ThreeBodyStorage_no2b.cc
    ThreeBodyStorage_mono.hh ThreeBodyStorage_mono.cc
    ThreeBodyStorage_iso.hh  ThreeBodyStorage_iso.cc
)
target_include_directories(
    ThreeBodyME
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
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
