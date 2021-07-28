add_library(cancer_biorobotscustmodules SHARED 
cancer_biorobots.cpp
)

target_link_libraries(cancer_biorobotscustmodules modules physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS cancer_biorobotscustmodules DESTINATION lib)
