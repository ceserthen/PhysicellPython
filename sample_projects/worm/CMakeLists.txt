add_subdirectory(custom_modules)

add_executable(worm main.cpp)

target_include_directories(worm PUBLIC
                          ${PROJECT_BINARY_DIR}
                          )
                          
target_link_libraries(worm PUBLIC wormcustmodules ${OpenMP_CXX_LIBRARIES})

install(DIRECTORY DESTINATION "worm/output")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/sample_projects/worm/config DESTINATION "worm")
install(TARGETS worm DESTINATION "worm")

IF(WIN32)
	install(TARGETS biofvm DESTINATION "worm")
	install(TARGETS physicellcore  DESTINATION "worm")
	install(TARGETS modules  DESTINATION "worm")
	install(TARGETS wormcustmodules DESTINATION "worm")
ENDIF()