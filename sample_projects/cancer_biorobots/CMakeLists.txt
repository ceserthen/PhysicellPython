add_subdirectory(custom_modules)

add_executable(cancer_biorobots main-cancer_biorobots.cpp)

target_include_directories(cancer_biorobots PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(cancer_biorobots PUBLIC cancer_biorobotscustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "cancer_biorobots/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/cancer_biorobots/config DESTINATION "cancer_biorobots")
install(TARGETS cancer_biorobots DESTINATION "cancer_biorobots")

IF(WIN32)
	install(TARGETS biofvm DESTINATION "cancer_biorobots")
	install(TARGETS physicellcore  DESTINATION "cancer_biorobots")
	install(TARGETS modules  DESTINATION "cancer_biorobots")
	install(TARGETS cancer_biorobotscustmodules DESTINATION "cancer_biorobots")
ENDIF()