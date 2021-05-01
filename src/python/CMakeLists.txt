cmake_minimum_required(VERSION 2.8)


add_library(pylib SHARED
     test_py.cpp)

target_link_libraries(pylib ${OpenMP_CXX_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_dependencies(pylib Boost)

# Tweaks the name of the library to match what Python expects
set_target_properties(pylib PROPERTIES SUFFIX .so)
set_target_properties(pylib PROPERTIES PREFIX "")

install(TARGETS pylib DESTINATION lib)
