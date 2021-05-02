cmake_minimum_required(VERSION 3.4)

pybind11_add_module(biofvmpy biofvm_py.cpp)

target_link_libraries(biofvmpy PRIVATE biofvm)

install(TARGETS biofvmpy DESTINATION lib)



