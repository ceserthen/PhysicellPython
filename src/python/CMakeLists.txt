cmake_minimum_required(VERSION 3.4)

add_subdirectory(BioFVM)

pybind11_add_module(physicell physicell_py.cpp)

target_link_libraries(physicell PRIVATE biofvmpy biofvm modules physicellcore )

install(TARGETS physicell DESTINATION lib)
