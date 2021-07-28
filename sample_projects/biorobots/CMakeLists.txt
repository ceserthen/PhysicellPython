add_subdirectory(custom_modules)

add_executable(biorobots main-biorobots.cpp)

target_include_directories(biorobots PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(biorobots PUBLIC biorobotscustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "biorobots/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/biorobots/config DESTINATION "biorobots")
install(TARGETS biorobots DESTINATION "biorobots")


IF(WIN32)
	install(TARGETS biofvm DESTINATION "biorobots")
	install(TARGETS physicellcore  DESTINATION "biorobots")
	install(TARGETS modules  DESTINATION "biorobots")
	install(TARGETS biorobotscustmodules DESTINATION "biorobots")
ENDIF()