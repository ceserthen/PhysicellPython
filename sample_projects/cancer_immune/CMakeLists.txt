cmake_minimum_required(VERSION 2.8)

add_subdirectory(custom_modules)

add_executable(heterogeneity main-heterogeneity.cpp)

target_include_directories(heterogeneity PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(heterogeneity PUBLIC hetcustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "heterogeneity/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/heterogeneity/config DESTINATION "heterogeneity")
install(TARGETS heterogeneity DESTINATION "heterogeneity")

IF(WIN32)
	install(TARGETS biofvm DESTINATION "cancer_immune")
	install(TARGETS physicellcore  DESTINATION "cancer_immune")
	install(TARGETS modules  DESTINATION "cancer_immune")
	install(TARGETS cancer_immunecustmodules DESTINATION "cancer_immune")
ENDIF()