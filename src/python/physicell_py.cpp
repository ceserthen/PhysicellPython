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

#include "../../pybind11/include/pybind11/pybind11.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>




#include <iostream>
#include <vector> 

#include "BioFVM/biofvm_py.h"
#include "../physicellcore/PhysiCell_cell.h"

namespace py = pybind11;



PYBIND11_MAKE_OPAQUE(std::vector<std::string, std::allocator<std::string>>);
using StringList = std::vector<std::string, std::allocator<std::string>>;

PYBIND11_MAKE_OPAQUE(std::vector<double, std::allocator<double>>);
using DoubleList = std::vector<double, std::allocator<double>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);


PYBIND11_MODULE(physicell, p) 
{
//BioFVM Modules
    py::module_ m  = p.def_submodule("BioFVM", "BioFVM submodule of Physicell"); 
    
    py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::bind_vector<std::vector<std::string>>(m, "VectorString");
    
    //TODO: Not currently working DAB 05-03-21
    //py::bind_vector<std::vector<std::vector<double>>>(m, "VectorVectorDouble");
    
    
//Microenvironment class
    py::class_<BioFVM::Microenvironment>(m, "Microenvironmentleg")
        //constructors
        .def(py::init<std::string>())
        .def(py::init<>())
        
        //attributes
        .def_readwrite("name", &BioFVM::Microenvironment::name)
        .def_readwrite("spatial_units", &BioFVM::Microenvironment::spatial_units)
        .def_readwrite("time_units", &BioFVM::Microenvironment::time_units)
        
        //diffusing entities
        .def_readwrite("density_names", &BioFVM::Microenvironment::density_names) 
        
        //coefficients
        .def_readwrite("diffusion_coefficients", &BioFVM::Microenvironment::diffusion_coefficients)
        .def_readwrite("decay_rates", &BioFVM::Microenvironment::decay_rates)
        .def_readwrite("supply_target_densities_times_supply_rates", &BioFVM::Microenvironment::supply_target_densities_times_supply_rates)
        .def_readwrite("supply_rates", &BioFVM::Microenvironment::supply_rates)
        .def_readwrite("uptake_rates", &BioFVM::Microenvironment::uptake_rates)
        
        .def("update_rates", static_cast<void (BioFVM::Microenvironment::*)(void)>(&BioFVM::Microenvironment::update_rates), "Update the rate coefficients")
        
        //TODO Function pointers for solvers
        //     void (*diffusion_decay_solver)( Microenvironment&, double); 
        // void (*bulk_supply_rate_function)( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination );
        // void (*bulk_supply_target_densities_function)( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination );
        // void (*bulk_uptake_rate_function)( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination );
        
        /*! functions to simplify size queries */ 
        .def("number_of_densities", static_cast<unsigned int (BioFVM::Microenvironment::*)(void)>(&BioFVM::Microenvironment::number_of_densities), "number_of_densities")
        .def("number_of_voxels", static_cast<unsigned int (BioFVM::Microenvironment::*)(void)>(&BioFVM::Microenvironment::number_of_voxels), "number_of_voxels")
        .def("number_of_voxel_faces", static_cast<unsigned int (BioFVM::Microenvironment::*)(void)>(&BioFVM::Microenvironment::number_of_voxel_faces), "number_of_voxel_faces")
        
        .def("auto_choose_diffusion_decay_solver", static_cast<void (BioFVM::Microenvironment::*)(void)>(&BioFVM::Microenvironment::auto_choose_diffusion_decay_solver), "auto_choose_diffusion_decay_solver")
        
        .def("resize_voxels", static_cast<void (BioFVM::Microenvironment::*)(int)>(&BioFVM::Microenvironment::resize_voxels), "resize_voxels: enter the new number of voxels", py::arg("new_number_of_voxels"))
        
        .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(int, int, int)>(&BioFVM::Microenvironment::resize_space), "resize_space", py::arg("x_nodes"), py::arg("y_nodes"), py::arg("z_nodes"))
        .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(double,double,double,double,double,double,int,int,int)> (&BioFVM::Microenvironment::resize_space), "resize_space", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("x_nodes"), py::arg("y_nodes"), py::arg("z_nodes"))
        .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(double,double,double,double,double,double,double,double,double)> (&BioFVM::Microenvironment::resize_space), "resize_space", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("dx_new"), py::arg("dy_new"), py::arg("dz_new"))
        
        .def("resize_space_uniform", static_cast<void (BioFVM::Microenvironment::*)(double,double,double,double,double,double,double)> (&BioFVM::Microenvironment::resize_space_uniform), "resize_space_uniform", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("dx_new"))
        
        //simulation methods
        .def("simulate_diffusion_decay", static_cast<void (BioFVM::Microenvironment::*)(double)>(&BioFVM::Microenvironment::simulate_diffusion_decay), "advance the diffusion-decay solver by dt time", py::arg("dt"))
        
        ;
   
        
        
        
//Microenvironment_py Class
    py::class_<BioFVM_py::Microenvironment_py, BioFVM::Microenvironment>(m, "Microenvironment")
    //constructors    
        .def(py::init<>())
        .def(py::init<BioFVM::Microenvironment_Options>(), py::arg("Microenvironment_Options"))
        .def(py::init<std::string>(), py::arg("xml_config_file"))
        .def(py::init<pugi::xml_node>(), py::arg("xml_node"))
    
    //XML readers
        
        .def("setup_microenvironment_from_file", static_cast<void (BioFVM_py::Microenvironment_py::*)(std::string)>(&BioFVM_py::Microenvironment_py::setup_microenvironment_from_file), "Load microenvironment settings from an XML file ", py::arg("xml_config_file"))
        
        .def("setup_microenvironment_from_XML_node", static_cast<bool (BioFVM_py::Microenvironment_py::*)(pugi::xml_node)>(&BioFVM_py::Microenvironment_py::setup_microenvironment_from_XML_node), "Load microenvironment settings from an XML node ", py::arg("xml_node"))
        
        .def("initialize_microenvironment", static_cast<void (BioFVM_py::Microenvironment_py::*)( void)>(&BioFVM_py::Microenvironment_py::initialize_microenvironment), "initialize microenvironment from preloaded xml configuration")
        
        .def("display_info", static_cast<void (BioFVM_py::Microenvironment_py::*)( void)>(&BioFVM_py::Microenvironment_py::display_info), "display information about the microenvironment")
        
        .def("resize_densities", static_cast<void (BioFVM_py::Microenvironment_py::*)( int)>(&BioFVM_py::Microenvironment_py::resize_densities), "resize_densities")
        
        .def("add_density", static_cast<void (BioFVM_py::Microenvironment_py::*)( void)>(&BioFVM_py::Microenvironment_py::add_density), "add_density")
        .def("add_density", static_cast<void (BioFVM_py::Microenvironment_py::*)( std::string, std::string)>(&BioFVM_py::Microenvironment_py::add_density), "add_density", py::arg("name"), py::arg("units"))
        .def("add_density", static_cast<void (BioFVM_py::Microenvironment_py::*)( std::string, std::string, double, double)>(&BioFVM_py::Microenvironment_py::add_density), "add_density", py::arg("name"), py::arg("units"), py::arg("diffusion_constant"), py::arg("decay_rate"))
        
        .def("set_density", static_cast<void (BioFVM_py::Microenvironment_py::*)( int, std::string, std::string)>(&BioFVM_py::Microenvironment_py::set_density), "set_density", py::arg("index"), py::arg("name"), py::arg("units"))
        .def("set_density", static_cast<void (BioFVM_py::Microenvironment_py::*)( int, std::string, std::string, double, double)>(&BioFVM_py::Microenvironment_py::set_density), "set_density", py::arg("index"), py::arg("name"), py::arg("units"), py::arg("diffusion_constant"), py::arg("decay_rate"))
        
        ;
        
        
//Core Modules
        py::module_ pcore  = p.def_submodule("core", "core modules for physicell"); 
        
//Cell Container
        py::class_<PhysiCell::Cell_Container>(pcore, "Cell_Container")
            //contructors
            .def(py::init<>())
            
            //attributes
            .def_readwrite("num_divisions_in_current_step", &PhysiCell::Cell_Container::num_divisions_in_current_step)
            .def_readwrite("num_deaths_in_current_step", &PhysiCell::Cell_Container::num_deaths_in_current_step)
            
            .def_readwrite("last_diffusion_time", &PhysiCell::Cell_Container::last_diffusion_time)
            .def_readwrite("last_cell_cycle_time", &PhysiCell::Cell_Container::last_cell_cycle_time)
            .def_readwrite("last_mechanics_time", &PhysiCell::Cell_Container::last_mechanics_time)
            
            //operators
            .def("initialize", static_cast<void (PhysiCell::Cell_Container::*) (double, double, double, double, double, double, double)>(&PhysiCell::Cell_Container::initialize), "Initialize the cell container by providing the ranges", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("voxel_size"))
            
            .def("initialize", static_cast<void (PhysiCell::Cell_Container::*) (double, double, double, double, double, double, double, double, double)>(&PhysiCell::Cell_Container::initialize), "Initialize the cell container by providing the ranges", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("dx"), py::arg("dy"), py::arg("dz"))
            
            .def("update_all_cells", static_cast<void (PhysiCell::Cell_Container::*) (double)>(&PhysiCell::Cell_Container::update_all_cells), "Update all cells", py::arg("t"))
            .def("update_all_cells", static_cast<void (PhysiCell::Cell_Container::*) (double, double)>(&PhysiCell::Cell_Container::update_all_cells), "Update all cells", py::arg("t"), py::arg("dt"))
            .def("update_all_cells", static_cast<void (PhysiCell::Cell_Container::*) (double, double, double)>(&PhysiCell::Cell_Container::update_all_cells), "Update all cells", py::arg("t"), py::arg("phenotype_dt"), py::arg("mechanics_dt"))
            .def("update_all_cells", static_cast<void (PhysiCell::Cell_Container::*) (double, double, double, double)>(&PhysiCell::Cell_Container::update_all_cells), "Update all cells", py::arg("t"), py::arg("phenotype_dt"), py::arg("mechanics_dt"), py::arg("diffusion_dt"))
            
            .def("register_agent", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*)> (&PhysiCell::Cell_Container::register_agent), "register an agent", py::arg("agent"))
            
            .def("add_agent_to_outer_voxel", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*)> (&PhysiCell::Cell_Container::add_agent_to_outer_voxel), "add an agent to an outer voxel", py::arg("agent"))
            
            .def("remove_agent", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*)> (&PhysiCell::Cell_Container::remove_agent), "remove an agent", py::arg("agent"))
            
            .def("remove_agent_from_voxel", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*, int)> (&PhysiCell::Cell_Container::remove_agent_from_voxel), "remove an agent from a voxel", py::arg("agent"), py::arg("voxel_index"))
        
            .def("add_agent_to_voxel", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*, int)> (&PhysiCell::Cell_Container::add_agent_to_voxel), "add an agent to a voxel", py::arg("agent"), py::arg("voxel_index"))
            
            .def("flag_cell_for_division", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*)> (&PhysiCell::Cell_Container::flag_cell_for_division), "flag a cell agent for division", py::arg("agent"))
            
            .def("flag_cell_for_removal", static_cast<void (PhysiCell::Cell_Container::*) (PhysiCell::Cell*)> (&PhysiCell::Cell_Container::flag_cell_for_removal), "flag a cell agent for removal", py::arg("agent"))
            
            .def("contain_any_cell", static_cast<bool (PhysiCell::Cell_Container::*) (int)> (&PhysiCell::Cell_Container::contain_any_cell), "contain_any_cell", py::arg("voxel_index"))
            ;

        //function for creating the cell container given a microenvironment
        pcore.def("create_cell_container_for_microenvironment", &PhysiCell::create_cell_container_for_microenvironment, py::arg("Microenvironment"), py::arg("mechanics_voxel_size"));

        

       
//Cell_Parameters
        py::class_<PhysiCell::Cell_Parameters>(pcore, "Cell_Parameters")
        //contructors
        .def(py::init<>())
        //attributes
        
        ;
        
        
//Cell_Definition
        py::class_<PhysiCell::Cell_Definition>(pcore, "Cell_Definition")
        //contructors
        .def(py::init<>())
        //attributes
        
        ;
        
        
//Cell_State
        py::class_<PhysiCell::Cell_State>(pcore, "Cell_State")
        //contructors
        .def(py::init<>())
        //attributes
        
        ;
        
//Cell  
        py::class_<PhysiCell::Cell>(pcore, "Cell")
        //contructors
        .def(py::init<>())
        //attributes
        
        ;
};
