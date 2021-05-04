cmake_minimum_required(VERSION 3.4)

add_library(physicellpy SHARED
     physicellcore_py.cpp)

target_link_libraries(physicellpy biofvmpy biofvm modules physicellcore )

install(TARGETS physicellpy DESTINATION lib)



