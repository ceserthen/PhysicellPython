add_library(hetcustmodules SHARED 
heterogeneity.cpp
)

target_link_libraries(hetcustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS hetcustmodules DESTINATION lib)
