cmake_minimum_required(VERSION 3.4)

pybind11_add_module(example biofvm_py.cpp)

install(TARGETS example DESTINATION lib)



