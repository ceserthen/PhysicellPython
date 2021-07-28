add_library(templatecustmodules SHARED 
custom.cpp
)

target_link_libraries(templatecustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS templatecustmodules DESTINATION lib)
