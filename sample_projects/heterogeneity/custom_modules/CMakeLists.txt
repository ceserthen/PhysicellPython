cmake_minimum_required(VERSION 2.8)

add_library(hetcustmodules SHARED 
heterogeneity.cpp
)

target_link_libraries(hetcustmodules modules core biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS hetcustmodules DESTINATION lib)
