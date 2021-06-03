cmake_minimum_required(VERSION 2.8)

add_library(physicellcore SHARED
PhysiCell_basic_signaling.cpp
PhysiCell_cell.cpp
PhysiCell_cell_container.cpp
PhysiCell_constants.cpp
PhysiCell_custom.cpp
#PhysiCell_digital_cell_line.cpp
PhysiCell_phenotype.cpp
PhysiCell_standard_models.cpp
PhysiCell_utilities.cpp
PhysiCell_settings.cpp
PhysiCell_pugixml.cpp
)

target_link_libraries(physicellcore biofvm ${OpenMP_CXX_LIBRARIES})

install(TARGETS physicellcore DESTINATION lib)

