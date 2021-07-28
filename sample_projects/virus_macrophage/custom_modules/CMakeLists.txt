add_library(virus_macrophagecustmodules SHARED 
custom.cpp
)

target_link_libraries(virus_macrophagecustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS virus_macrophagecustmodules DESTINATION lib)
