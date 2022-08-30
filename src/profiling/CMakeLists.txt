add_executable(comm232 comm232.cc)
target_link_libraries(comm232 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)

add_executable(comm223 comm223.cc)
target_link_libraries(comm223 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)
