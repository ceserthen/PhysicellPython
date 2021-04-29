cmake_minimum_required(VERSION 2.8)

add_library(wormcustmodules SHARED 
custom.cpp
)

target_link_libraries(wormcustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS wormcustmodules DESTINATION lib)
