cmake_minimum_required(VERSION 2.8)

add_library(templatecustmodules SHARED 
custom.cpp
)

target_link_libraries(templatecustmodules modules core biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS templatecustmodules DESTINATION lib)
