cmake_minimum_required(VERSION 2.8)

add_subdirectory(custom_modules)

add_executable(virus_macrophage main.cpp)

target_include_directories(virus_macrophage PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(virus_macrophage PUBLIC virus_macrophagecustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "virus_macrophage/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/virus_macrophage/config DESTINATION "virus_macrophage")
install(TARGETS virus_macrophage DESTINATION "virus_macrophage")

