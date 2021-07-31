/*
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#ifndef __physicellcore_py_h__
#define __physicellcore_py_h__

#include "../../../src/physicellcore/PhysiCell.h"
#include "../BioFVM/biofvm_py.h"
//TODO:Pathology
//TODO:Various Outputs
namespace PhysiCell{
	// pybind11 functions originally created by @rheiland for https://github.com/rheiland/PyPhysiCell
void create_cell_container2( double mechanics_voxel_size );
// pybind11
void update_all_cells(double t);
// pybind11
int get_num_cells( void );
// pybind11
std::vector<double> get_cells_pos2D( void );
};
namespace PhysiCellCore_py{

class Cell_Container_py : public PhysiCell::Cell_Container {
public:
    std::vector<PhysiCell::Cell*> *all_cells; 
    
    Cell_Container_py();
    Cell_Container_py(BioFVM::Microenvironment& m , double mechanics_voxel_size);
    
    
    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size);
    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz);
    
    void create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size ); 
    
    void update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ ); 
    void update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ );
    
    // making the create and save cell functions part of the 
    PhysiCell::Cell* create_cell( void );
    PhysiCell::Cell* create_cell( PhysiCell::Cell_Definition& cd );
    //need to add a die and divide function to compliment performance
    void delete_cell( int index );
    void delete_cell(PhysiCell::Cell* pDelete); 
    
    void kill_cell(PhysiCell::Cell* pDelete);
    PhysiCell::Cell* divide_cell(PhysiCell::Cell* pDivide);
    
    void save_all_cells_to_matlab( std::string filename ); //TODO
    void delete_cell_original( int index ); 
    
};

void add_PhysiCell_cells_to_open_xml_pugi_py( pugi::xml_document& xml_dom, std::string filename_base, Microenvironment& M, Cell_Container_py& CellCon  );

void save_PhysiCell_to_MultiCellDS_xml_pugi_py( std::string filename_base , BioFVM_py::Microenvironment_py& M , Cell_Container_py& CellCon, double current_simulation_time);


    
}//end namespace PhysiCellCore_py


#endif
