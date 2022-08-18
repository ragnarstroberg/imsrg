add_executable(comm232 comm232.cc)
target_link_libraries(comm232 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)
