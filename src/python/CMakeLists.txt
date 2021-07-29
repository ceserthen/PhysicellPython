add_subdirectory(BioFVM)
add_subdirectory(physicellcore)

pybind11_add_module(physicell physicell_py.cpp)

target_link_libraries(physicell PRIVATE physicellpy biofvmpy biofvm modules physicellcore )

install(TARGETS physicell DESTINATION lib)
