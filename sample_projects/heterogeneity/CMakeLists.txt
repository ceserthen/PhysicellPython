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
	install(TARGETS biofvm DESTINATION "heterogeneity")
	install(TARGETS physicellcore  DESTINATION "heterogeneity")
	install(TARGETS modules  DESTINATION "heterogeneity")
	install(TARGETS hetcustmodules DESTINATION "heterogeneity")
ENDIF()