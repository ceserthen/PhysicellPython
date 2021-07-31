cmake_minimum_required(VERSION 3.4)

add_library(biorobotspy SHARED
     biorobots_py.cpp)

target_link_libraries(biorobotspy physicellpy biofvmpy biofvm modules physicellcore )

install(TARGETS biorobotspy DESTINATION lib)



