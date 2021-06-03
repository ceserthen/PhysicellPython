cmake_minimum_required(VERSION 2.8)

add_library(biorobotscustmodules SHARED 
biorobots.cpp
)

target_link_libraries(biorobotscustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS biorobotscustmodules DESTINATION lib)
