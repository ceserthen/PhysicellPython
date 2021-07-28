cmake_minimum_required(VERSION 3.4)

add_library(biofvmpy SHARED
     biofvm_py.cpp)

target_link_libraries(biofvmpy biofvm modules physicellcore )

install(TARGETS biofvmpy DESTINATION lib)



