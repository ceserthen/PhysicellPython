cmake_minimum_required(VERSION 2.8)

add_subdirectory(custom_modules)

add_executable(cancer_immune main-cancer_immune_3D.cpp)

target_include_directories(cancer_immune PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(cancer_immune PUBLIC cancer_immunecustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "cancer_immune/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/cancer_immune/config DESTINATION "cancer_immune")
install(TARGETS cancer_immune DESTINATION "cancer_immune")

