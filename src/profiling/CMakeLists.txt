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

add_executable(comm133 comm133.cc)
target_link_libraries(comm133 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)

add_executable(comm222pphh comm222pphh.cc)
target_link_libraries(comm222pphh 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)

add_executable(comm222ph comm222ph.cc)
target_link_libraries(comm222ph 
    IMSRG
    Parameters
    PhysicalConstants
    Version
)
