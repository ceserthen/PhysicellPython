cmake_minimum_required(VERSION 2.8)

add_library(virus_macrophagecustmodules SHARED 
custom.cpp
)

target_link_libraries(virus_macrophagecustmodules modules core biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS virus_macrophagecustmodules DESTINATION lib)
