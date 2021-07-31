add_subdirectory(BioFVM)
add_subdirectory(physicellcore)
add_subdirectory(tests)

pybind11_add_module(physicell physicell_py.cpp)

target_link_libraries(physicell PRIVATE biorobotspy physicellpy biofvmpy biofvm modules physicellcore )

install(TARGETS physicell DESTINATION lib)
