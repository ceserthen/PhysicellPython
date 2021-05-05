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
#include "../physicellcore/PhysiCell.h"

namespace py = pybind11;



PYBIND11_MAKE_OPAQUE(std::vector<std::string, std::allocator<std::string>>);
using StringList = std::vector<std::string, std::allocator<std::string>>;

PYBIND11_MAKE_OPAQUE(std::vector<double, std::allocator<double>>);
using DoubleList = std::vector<double, std::allocator<double>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);


PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Cell*, std::allocator<PhysiCell::Cell*>>);
using CellList = std::vector<PhysiCell::Cell*, std::allocator<PhysiCell::Cell*>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PhysiCell::Cell*>>);


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
        .def_readwrite("o2_hypoxic_threshold", &PhysiCell::Cell_Parameters::o2_hypoxic_threshold, "value at which hypoxic signaling starts")
        .def_readwrite("o2_hypoxic_response", &PhysiCell::Cell_Parameters::o2_hypoxic_response, "value at which omics changes are observed")
        .def_readwrite("o2_hypoxic_saturation", &PhysiCell::Cell_Parameters::o2_hypoxic_saturation, "value at which hypoxic signalign saturates (o2_hypoxic_saturation < o2_hypoxic_threshold)")
        
        .def_readwrite("o2_proliferation_saturation", &PhysiCell::Cell_Parameters::o2_proliferation_saturation, "value at which extra o2 does not increase proliferation")
        .def_readwrite("o2_proliferation_threshold", &PhysiCell::Cell_Parameters::o2_proliferation_threshold, "value at which o2 is sufficient for proliferation")
        
        .def_readwrite("o2_reference", &PhysiCell::Cell_Parameters::o2_reference, "physioxic reference value, in the linked reference Phenotype (o2_proliferation_threshold < o2_reference < o2_proliferation_saturation)")
        
        .def_readwrite("o2_necrosis_threshold", &PhysiCell::Cell_Parameters::o2_necrosis_threshold, "value at which cells start experiencing necrotic death")
        .def_readwrite("o2_necrosis_max", &PhysiCell::Cell_Parameters::o2_necrosis_max, "value at which necrosis reaches its maximum rate (o2_necrosis_max < o2_necrosis_threshold)")
        
        .def_readwrite("pReference_live_phenotype", &PhysiCell::Cell_Parameters::pReference_live_phenotype)
        
        .def_readwrite("max_necrosis_rate", &PhysiCell::Cell_Parameters::max_necrosis_rate)
        .def_readwrite("necrosis_type", &PhysiCell::Cell_Parameters::necrosis_type)
        ;
        
        
//Cell_Definition
        py::class_<PhysiCell::Cell_Definition>(pcore, "Cell_Definition")
        //contructors
        .def(py::init<>())
        //.def(py::init<PhysiCell::Cell_Definition&)>()
        
        //operatos
        //.def(py::self = py::self)
        
        //attributes
        .def_readwrite("type", &PhysiCell::Cell_Definition::type)
        .def_readwrite("name", &PhysiCell::Cell_Definition::name)
        
        .def_readwrite("pMicroenvironment", &PhysiCell::Cell_Definition::pMicroenvironment)
        
        .def_readwrite("parameters", &PhysiCell::Cell_Definition::parameters)
        .def_readwrite("custom_data", &PhysiCell::Cell_Definition::custom_data)
        .def_readwrite("functions", &PhysiCell::Cell_Definition::functions)
        .def_readwrite("phenotype", &PhysiCell::Cell_Definition::phenotype)
        ;
        
        
//Cell_State
        py::class_<PhysiCell::Cell_State>(pcore, "Cell_State")
        //contructors
        .def(py::init<>())
        //attributes
        .def_readwrite("attached_cells", &PhysiCell::Cell_State::attached_cells)
        
        .def_readwrite("neighbors", &PhysiCell::Cell_State::neighbors)
        .def_readwrite("orientation", &PhysiCell::Cell_State::orientation)
        
        .def_readwrite("simple_pressure", &PhysiCell::Cell_State::simple_pressure)
        
        .def("number_of_attached_cells", static_cast<int (PhysiCell::Cell_State::*) (void)> (&PhysiCell::Cell_State::number_of_attached_cells), "Return the number of attached cells")
        ;
        
//Cell  
        py::class_<PhysiCell::Cell>(pcore, "Cell")
        //contructors
        .def(py::init<>())
        //attributes
        .def_readwrite("type_name", &PhysiCell::Cell::type_name)
        
        .def_readwrite("custom_data", &PhysiCell::Cell::custom_data)
        .def_readwrite("parameters", &PhysiCell::Cell::parameters)
        .def_readwrite("functions", &PhysiCell::Cell::functions)
        
        .def_readwrite("state", &PhysiCell::Cell::state)
        .def_readwrite("phenotype", &PhysiCell::Cell::phenotype)
        
        .def_readwrite("is_out_of_domain", &PhysiCell::Cell::is_out_of_domain)
        .def_readwrite("is_movable", &PhysiCell::Cell::is_movable)
        
        .def_readwrite("displacement", &PhysiCell::Cell::displacement)
        
        //operators
        .def("update_motility_vector", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::update_motility_vector), "update the motility vector", py::arg("dt"))
        
        .def("advance_bundled_phenotype_functions", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::advance_bundled_phenotype_functions), "advance bundled phenotype functions", py::arg("dt"))
        
        
        .def("add_potentials", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::add_potentials), "add potentails")
        
        .def("set_previous_velocity", static_cast<void (PhysiCell::Cell::*) (double, double, double)> (&PhysiCell::Cell::set_previous_velocity), "set_previous_velocity", py::arg("xV"), py::arg("yV"), py::arg("zV"))
        
        .def("get_current_mechanics_voxel_index", static_cast<int (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::get_current_mechanics_voxel_index), "get_current_mechanics_voxel_index")
        
        .def("turn_off_reactions", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::turn_off_reactions), "turn_off_reactions")
        
        
        .def("flag_for_division", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::flag_for_division), "flag_for_division")
        
        .def("flag_for_removal", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::flag_for_removal), "flag_for_removal")
        
        
        .def("start_death", static_cast<void (PhysiCell::Cell::*) (int)> (&PhysiCell::Cell::start_death), "start death given a ceartain death model", py::arg("death_model_index"))
        
        .def("lyse_cell", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::lyse_cell))
        
        
        .def("divide", static_cast<PhysiCell::Cell* (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::divide), "divide the cell into a new cell")
        
        .def("die", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::die), "the cells dies")
        
        .def("step", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::step), "step the cell", py::arg("dt"))
        
        
        .def("assign_position", static_cast<bool (PhysiCell::Cell::*) (std::vector<double>)> (&PhysiCell::Cell::assign_position), "assign_position")
        
        .def("assign_position", static_cast<bool (PhysiCell::Cell::*) (double, double, double)> (&PhysiCell::Cell::assign_position), "assign_position")
        
        .def("set_total_volume", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::set_total_volume), "set_total_volume")
        
        
        .def("get_total_volume", static_cast<double& (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::get_total_volume), "set_total_volume")
        
        
        .def("set_target_volume", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::set_target_volume), "set_target_volume")
        
        .def("set_target_radius", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::set_target_radius), "set_target_radius")
        
        .def("set_radius", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::set_radius), "set_radius")
        
        
        .def("update_position", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::update_position), "update_position", py::arg("dt"))
        
        
        .def("copy_function_pointers", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::copy_function_pointers), "copy_function_pointers")
        
        
        .def("assign_orientation", static_cast<void (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::assign_orientation), "assign_orientation")
        
        
        .def("update_voxel_in_container", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::update_voxel_in_container), "update_voxel_in_container")
        
        .def("copy_data", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::copy_data), "copy_data")
        
        
        .def("ingest_cell", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::ingest_cell), "ingest_cell")
        
        .def("attach_cell", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::attach_cell), "attach_cell")
        
        .def("detach_cell", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Cell*)> (&PhysiCell::Cell::detach_cell), "detach_cell")
        
        .def("remove_all_attached_cells", static_cast<void (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::remove_all_attached_cells), "remove_all_attached_cells")
        
        
        .def("set_phenotype", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Phenotype&)> (&PhysiCell::Cell::set_phenotype), "set_phenotype")
        
        .def("update_radius", static_cast<void (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::update_radius), "update_radius")
        
        .def("get_container", static_cast<PhysiCell::Cell_Container* (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::get_container), "get_container")
        
        .def("cells_in_my_container", static_cast<std::vector<PhysiCell::Cell*>& (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::cells_in_my_container), "cells_in_my_container")
        
        .def("nearby_cells", static_cast<std::vector<PhysiCell::Cell*> (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::nearby_cells), "nearby_cells")
        
        .def("nearby_interacting_cells", static_cast<std::vector<PhysiCell::Cell*> (PhysiCell::Cell::*) ( void)> (&PhysiCell::Cell::nearby_interacting_cells), "nearby_interacting_cells")
        
        .def("convert_to_cell_definition", static_cast<void (PhysiCell::Cell::*) ( PhysiCell::Cell_Definition&)> (&PhysiCell::Cell::convert_to_cell_definition), "convert_to_cell_definition")
        
        ;
};
