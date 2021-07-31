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
//#include <pybind11/functional.h>




#include <iostream>
#include <vector> 

#include "BioFVM/biofvm_py.h"
#include "../physicellcore/PhysiCell.h"
#include "physicellcore/physicellcore_py.h"

namespace py = pybind11;



PYBIND11_MAKE_OPAQUE(std::vector<std::string, std::allocator<std::string>>);
using StringList = std::vector<std::string, std::allocator<std::string>>;

PYBIND11_MAKE_OPAQUE(std::vector<unsigned int, std::allocator<unsigned int>>);
using UIntList = std::vector<unsigned int, std::allocator<unsigned int>>;

PYBIND11_MAKE_OPAQUE(std::vector<double, std::allocator<double>>);
using DoubleList = std::vector<double, std::allocator<double>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);


PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Cell*, std::allocator<PhysiCell::Cell*>>);
using CellPointerList = std::vector<PhysiCell::Cell*, std::allocator<PhysiCell::Cell*>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PhysiCell::Cell*>>);

PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Cell, std::allocator<PhysiCell::Cell>>);
using CellList = std::vector<PhysiCell::Cell, std::allocator<PhysiCell::Cell>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PhysiCell::Cell>>);


PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Cycle_Model*, std::allocator<PhysiCell::Cycle_Model*>>);
using CycleModelList = std::vector<PhysiCell::Cycle_Model*, std::allocator<PhysiCell::Cycle_Model*>>;

PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Death_Parameters*, std::allocator<PhysiCell::Death_Parameters*>>);
using DeathParametersList = std::vector<PhysiCell::Death_Parameters*, std::allocator<PhysiCell::Death_Parameters*>>;

PYBIND11_MAKE_OPAQUE(std::vector<PhysiCell::Phase, std::allocator<PhysiCell::Phase>>);
using PhaseList = std::vector<PhysiCell::Phase, std::allocator<PhysiCell::Phase>>;

PYBIND11_MAKE_OPAQUE(std::vector< std::vector<PhysiCell::Phase_Link>>);

PYBIND11_MODULE(physicell, p) 
{
//BioFVM Modules
    py::module_ m  = p.def_submodule("BioFVM", "BioFVM submodule of Physicell"); 
    
    py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::bind_vector<std::vector<std::string>>(m, "VectorString");
    
    //TODO: Not currently working DAB 05-03-21
    //py::bind_vector<std::vector<std::vector<double>>>(m, "VectorVectorDouble");
    
    py::class_<pugi::xml_document>(m, "xml_document_pugi")
//contructors
    .def(py::init<>())
    ;
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
        
        .def("resize_densities", static_cast<void (BioFVM::Microenvironment::*)(int)> (&BioFVM::Microenvironment::resize_densities), "resize_densities", py::arg("new_size"))
        .def("add_density", static_cast<void (BioFVM::Microenvironment::*)(void)> (&BioFVM::Microenvironment::add_density), "add_density")
        .def("add_density", static_cast<void (BioFVM::Microenvironment::*)(std::string, std::string)> (&BioFVM::Microenvironment::add_density), "add_density", py::arg("name"), py::arg("units"))
        .def("add_density", static_cast<void (BioFVM::Microenvironment::*)(std::string, std::string, double, double)> (&BioFVM::Microenvironment::add_density), "add_density", py::arg("name"), py::arg("units"), py::arg("diffusion_constan"), py::arg("decay_rate"))
        
        .def("set_density", static_cast<void (BioFVM::Microenvironment::*)(int, std::string, std::string)> (&BioFVM::Microenvironment::set_density), "set_density", py::arg("index"), py::arg("name"), py::arg("units"))
        .def("set_density", static_cast<void (BioFVM::Microenvironment::*)(int, std::string, std::string, double, double)> (&BioFVM::Microenvironment::set_density), "set_density", py::arg("index"),py::arg("name"), py::arg("units"), py::arg("diffusion_constan"), py::arg("decay_rate"))

        .def("find_density_index", static_cast<int (BioFVM::Microenvironment::*)(std::string)>(&BioFVM::Microenvironment::find_density_index), "add_density", py::arg("name"))
	
        .def("voxel_index", static_cast<int (BioFVM::Microenvironment::*)(int, int, int)>(&BioFVM::Microenvironment::voxel_index), "voxel_index", py::arg("i"), py::arg("j"), py::arg("k"))

        .def("cartesian_indices", static_cast<std::vector<unsigned int> (BioFVM::Microenvironment::*)(int)>(&BioFVM::Microenvironment::cartesian_indices), "cartesian_indices", py::arg("n"))
	
	// int nearest_voxel_index( std::vector<double>& position ); 
	// std::vector<unsigned int> nearest_cartesian_indices( std::vector<double>& position ); 
	// Voxel& nearest_voxel( std::vector<double>& position ); 
	// Voxel& voxels( int voxel_index );
	// std::vector<double>& nearest_density_vector( std::vector<double>& position );  
	// std::vector<double>& nearest_density_vector( int voxel_index );  

	// /*! access the density vector at  [ X(i),Y(j),Z(k) ] */
	// std::vector<double>& operator()( int i, int j, int k ); 
	// /*! access the density vector at  [ X(i),Y(j),0 ]  -- helpful for 2-D problems */
	// std::vector<double>& operator()( int i, int j );  
	// /*! access the density vector at [x,y,z](n) */
	// std::vector<double>& operator()( int n );  
	
	// std::vector<gradient>& gradient_vector(int i, int j, int k); 
	// std::vector<gradient>& gradient_vector(int i, int j ); 
	// std::vector<gradient>& gradient_vector(int n );  
	
	// std::vector<gradient>& nearest_gradient_vector( std::vector<double>& position ); 

	// void compute_all_gradient_vectors( void ); 
	// void compute_gradient_vector( int n );  
	// void reset_all_gradient_vectors( void ); 
	
	// /*! access the density vector at  [ X(i),Y(j),Z(k) ] */
	// std::vector<double>& density_vector( int i, int j, int k ); 
	// /*! access the density vector at  [ X(i),Y(j),0 ]  -- helpful for 2-D problems */
	// std::vector<double>& density_vector( int i, int j ); 
	// /*! access the density vector at [x,y,z](n) */
	// std::vector<double>& density_vector( int n ); 

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
        py::class_<PhysiCell::Cell_Container>(pcore, "Cell_Container_legacy")
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
            
            //.def("step", static_cast<void (PhysiCell::Cell::*) (double)> (&PhysiCell::Cell::step), "step the cell", py::arg("dt"))
            
            
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
            
            
            //.def("set_phenotype", static_cast<void (PhysiCell::Cell::*) (PhysiCell::Phenotype&)> (&PhysiCell::Cell::set_phenotype), "set_phenotype")
            
            //.def("update_radius", static_cast<void (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::update_radius), "update_radius")
            
            .def("get_container", static_cast<PhysiCell::Cell_Container* (PhysiCell::Cell::*) ()> (&PhysiCell::Cell::get_container), "get_container")
            
            .def("cells_in_my_container", static_cast<std::vector<PhysiCell::Cell*>& (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::cells_in_my_container), "cells_in_my_container")
            
            .def("nearby_cells", static_cast<std::vector<PhysiCell::Cell*> (PhysiCell::Cell::*) (void)> (&PhysiCell::Cell::nearby_cells), "nearby_cells")
            
            .def("nearby_interacting_cells", static_cast<std::vector<PhysiCell::Cell*> (PhysiCell::Cell::*) ( void)> (&PhysiCell::Cell::nearby_interacting_cells), "nearby_interacting_cells")
            
            .def("convert_to_cell_definition", static_cast<void (PhysiCell::Cell::*) ( PhysiCell::Cell_Definition&)> (&PhysiCell::Cell::convert_to_cell_definition), "convert_to_cell_definition")
            
        ;
        
        
    //Phase
    py::class_<PhysiCell::Phase>(pcore, "Phase")
        //constructor
        .def(py::init<>())
        
        //attributes
        .def_readwrite("index", &PhysiCell::Phase::index)
        .def_readwrite("code", &PhysiCell::Phase::code)
        .def_readwrite("name", &PhysiCell::Phase::name)
        .def_readwrite("division_at_phase_exit", &PhysiCell::Phase::division_at_phase_exit)
        .def_readwrite("removal_at_phase_exit", &PhysiCell::Phase::removal_at_phase_exit)
        //TODO: Entry Function
        
        ;
        
    //Phase_link
    py::class_<PhysiCell::Phase_Link>(pcore, "Phase_Link")
        //constructor
        .def(py::init<>())
        
        //attributes
        .def_readwrite("start_phase_index", &PhysiCell::Phase_Link::start_phase_index)    
        .def_readwrite("end_phase_index", &PhysiCell::Phase_Link::end_phase_index)  
        .def_readwrite("fixed_duration", &PhysiCell::Phase_Link::fixed_duration)  
        
        ;
        
    //Cycle_Data
    py::class_<PhysiCell::Cycle_Data>(pcore, "Cycle_Data")
        //constructor
        .def(py::init<>())
        
        //attributes
        .def_readwrite("pCycle_Model", &PhysiCell::Cycle_Data::pCycle_Model)
        .def_readwrite("time_units", &PhysiCell::Cycle_Data::time_units)
        .def_readwrite("transition_rates", &PhysiCell::Cycle_Data::transition_rates)
        .def_readwrite("current_phase_index", &PhysiCell::Cycle_Data::current_phase_index)
        .def_readwrite("elapsed_time_in_phase", &PhysiCell::Cycle_Data::elapsed_time_in_phase)
        
        //operators
        .def("current_phase", static_cast<PhysiCell::Phase& (PhysiCell::Cycle_Data::*) (void)> (&PhysiCell::Cycle_Data::current_phase), "current_phase")
        
        .def("sync_to_cycle_model", static_cast<void (PhysiCell::Cycle_Data::*) (void)> (&PhysiCell::Cycle_Data::sync_to_cycle_model), "sync_to_cycle_model")
        
        .def("transition_rate", static_cast<double& (PhysiCell::Cycle_Data::*) (int,int)> (&PhysiCell::Cycle_Data::transition_rate), "transition_rate", py::arg("start_phase_index"), py::arg("end_phase_index"))
        
        .def("exit_rate", static_cast<double& (PhysiCell::Cycle_Data::*) (int)> (&PhysiCell::Cycle_Data::exit_rate), "exit_rate")
        
        ;
        
        
    //Cycle_Model
    py::class_<PhysiCell::Cycle_Model>(pcore, "Cycle_Model")
        //constructor
        .def(py::init<>())
        
        //attributes
        .def_readwrite("name", &PhysiCell::Cycle_Model::name)
        .def_readwrite("code", &PhysiCell::Cycle_Model::code)
        .def_readwrite("phases", &PhysiCell::Cycle_Model::phases)
        .def_readwrite("phase_links", &PhysiCell::Cycle_Model::phase_links)
        .def_readwrite("default_phase_index", &PhysiCell::Cycle_Model::default_phase_index)
        .def_readwrite("data", &PhysiCell::Cycle_Model::data)
        
        //operators
        .def("advance_model", static_cast<void (PhysiCell::Cycle_Model::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Cycle_Model::advance_model), "advance_model")
        
        .def("add_phase", static_cast<int (PhysiCell::Cycle_Model::*) (int, std::string)> (&PhysiCell::Cycle_Model::add_phase), "add_phase", py::arg("start_index"), py::arg("end_index"))
        
        //TODO: put together the add hase linked
        
        .def("find_phase_index", static_cast<int (PhysiCell::Cycle_Model::*) (int)> (&PhysiCell::Cycle_Model::find_phase_index), "find_phase_index", py::arg("code"))
        .def("find_phase_index", static_cast<int (PhysiCell::Cycle_Model::*) (std::string)> (&PhysiCell::Cycle_Model::find_phase_index), "find_phase_index", py::arg("name"))
        
         .def("transition_rate", static_cast<double& (PhysiCell::Cycle_Model::*) (int, int)> (&PhysiCell::Cycle_Model::transition_rate), "transition_rate", py::arg("start_index"), py::arg("end_index"))
         
         
        .def("phase_link", static_cast<PhysiCell::Phase_Link& (PhysiCell::Cycle_Model::*) (int, int)> (&PhysiCell::Cycle_Model::phase_link), "phase_link", py::arg("start_index"), py::arg("end_index"))
        
        ;
        
    //Cycle
    py::class_<PhysiCell::Cycle>(pcore, "Cycle")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("pCycle_Model", &PhysiCell::Cycle::pCycle_Model)
        .def_readwrite("data", &PhysiCell::Cycle::data) 
        
        //operators
        .def("advance_cycle", static_cast<void (PhysiCell::Cycle::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Cycle::advance_cycle), "advance_cycle")
        
        .def("model", static_cast<PhysiCell::Cycle_Model& (PhysiCell::Cycle::*) (void)> (&PhysiCell::Cycle::model), "model")
        
        .def("current_phase", static_cast<PhysiCell::Phase& (PhysiCell::Cycle::*) (void)> (&PhysiCell::Cycle::current_phase), "current_phase")
        
        .def("current_phase_index", static_cast<int& (PhysiCell::Cycle::*) (void)> (&PhysiCell::Cycle::current_phase_index), "current_phase_index")
        
        
        .def("sync_to_cycle_model", static_cast<void (PhysiCell::Cycle::*) (PhysiCell::Cycle_Model&)> (&PhysiCell::Cycle::sync_to_cycle_model), "sync_to_cycle_model")
        
        
        ;
        
    //Death_Parameter
    py::class_<PhysiCell::Death_Parameters>(pcore, "Death_Parameters")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("time_units", &PhysiCell::Death_Parameters::time_units)
        .def_readwrite("unlysed_fluid_change_rate", &PhysiCell::Death_Parameters::unlysed_fluid_change_rate)
        .def_readwrite("lysed_fluid_change_rate", &PhysiCell::Death_Parameters::lysed_fluid_change_rate)
        .def_readwrite("cytoplasmic_biomass_change_rate", &PhysiCell::Death_Parameters::cytoplasmic_biomass_change_rate)
        .def_readwrite("nuclear_biomass_change_rate", &PhysiCell::Death_Parameters::nuclear_biomass_change_rate)
        .def_readwrite("calcification_rate", &PhysiCell::Death_Parameters::calcification_rate)
        .def_readwrite("relative_rupture_volume", &PhysiCell::Death_Parameters::relative_rupture_volume)
        ;
    //Death
    py::class_<PhysiCell::Death>(pcore, "Death")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("rates", &PhysiCell::Death::rates)
        .def_readwrite("models", &PhysiCell::Death::models)
        .def_readwrite("parameters", &PhysiCell::Death::parameters)
        
        .def_readwrite("dead", &PhysiCell::Death::dead)
        .def_readwrite("current_death_model_index", &PhysiCell::Death::current_death_model_index)
        
        //operators
        .def("add_death_model", static_cast<int (PhysiCell::Death::*) (double, PhysiCell::Cycle_Model*)> (&PhysiCell::Death::add_death_model), "add_death_model")
        
        .def("add_death_model", static_cast<int (PhysiCell::Death::*) (double, PhysiCell::Cycle_Model*, PhysiCell::Death_Parameters&)> (&PhysiCell::Death::add_death_model), "add_death_model")
        
        .def("find_death_model_index", static_cast<int (PhysiCell::Death::*) (int)> (&PhysiCell::Death::find_death_model_index), "find_death_model_index")
        .def("find_death_model_index", static_cast<int (PhysiCell::Death::*) (std::string)> (&PhysiCell::Death::find_death_model_index), "find_death_model_index")
        
        .def("check_for_death", static_cast<bool (PhysiCell::Death::*) (double)> (&PhysiCell::Death::check_for_death), "check_for_death")
        
        .def("trigger_death", static_cast<void (PhysiCell::Death::*) (int)> (&PhysiCell::Death::trigger_death), "trigger_death", py::arg("death_model_index"))
        
        
        .def("current_model", static_cast<PhysiCell::Cycle_Model& (PhysiCell::Death::*) (void)> (&PhysiCell::Death::current_model), "current_model")
        
        .def("current_parameters", static_cast<PhysiCell::Death_Parameters& (PhysiCell::Death::*) (void)> (&PhysiCell::Death::current_parameters), "current_parameters")
        
        ;
    //Volume
    py::class_<PhysiCell::Volume>(pcore, "Volume")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("total", &PhysiCell::Volume::total)
        .def_readwrite("solid", &PhysiCell::Volume::solid)
        .def_readwrite("fluid", &PhysiCell::Volume::fluid)
        .def_readwrite("fluid_fraction", &PhysiCell::Volume::fluid_fraction)
        
        .def_readwrite("nuclear", &PhysiCell::Volume::nuclear)
        .def_readwrite("nuclear_fluid", &PhysiCell::Volume::nuclear_fluid)
        .def_readwrite("nuclear_solid", &PhysiCell::Volume::nuclear_solid)
        
        .def_readwrite("cytoplasmic", &PhysiCell::Volume::cytoplasmic)
        .def_readwrite("cytoplasmic_fluid", &PhysiCell::Volume::cytoplasmic_fluid)
        .def_readwrite("cytoplasmic_solid", &PhysiCell::Volume::cytoplasmic_solid)
        
        .def_readwrite("calcified_fraction", &PhysiCell::Volume::calcified_fraction)
        
        .def_readwrite("cytoplasmic_to_nuclear_ratio", &PhysiCell::Volume::cytoplasmic_to_nuclear_ratio)
        
        .def_readwrite("rupture_volume", &PhysiCell::Volume::rupture_volume)
        
        .def_readwrite("cytoplasmic_biomass_change_rate", &PhysiCell::Volume::cytoplasmic_biomass_change_rate)
        .def_readwrite("nuclear_biomass_change_rate", &PhysiCell::Volume::nuclear_biomass_change_rate)
        .def_readwrite("fluid_change_rate", &PhysiCell::Volume::fluid_change_rate)
        
        .def_readwrite("calcification_rate", &PhysiCell::Volume::calcification_rate)
        
        .def_readwrite("target_solid_cytoplasmic", &PhysiCell::Volume::target_solid_cytoplasmic)
        .def_readwrite("target_solid_nuclear", &PhysiCell::Volume::target_solid_nuclear)
        .def_readwrite("target_fluid_fraction", &PhysiCell::Volume::target_fluid_fraction)
        
        .def_readwrite("target_cytoplasmic_to_nuclear_ratio", &PhysiCell::Volume::target_cytoplasmic_to_nuclear_ratio)
        
        .def_readwrite("relative_rupture_volume", &PhysiCell::Volume::relative_rupture_volume)
        
        //operators
        .def("divide", static_cast<void (PhysiCell::Volume::*) (void)> (&PhysiCell::Volume::divide), "divide")
        
        .def("multiply_by_ratio", static_cast<void (PhysiCell::Volume::*) (double)> (&PhysiCell::Volume::multiply_by_ratio), "multiply_by_ratio")
        ;     

    //Geometry
    py::class_<PhysiCell::Geometry>(pcore, "Geometry")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("radius", &PhysiCell::Geometry::radius)
        .def_readwrite("nuclear_radius", &PhysiCell::Geometry::nuclear_radius)
        .def_readwrite("surface_area", &PhysiCell::Geometry::surface_area)
        
        .def_readwrite("polarity", &PhysiCell::Geometry::polarity)
        
        //operators
        .def("update_radius", static_cast<void (PhysiCell::Geometry::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Geometry::update_radius), "update_radius")
        .def("update_nuclear_radius", static_cast<void (PhysiCell::Geometry::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Geometry::update_nuclear_radius), "update_nuclear_radius")
        .def("update_surface_area", static_cast<void (PhysiCell::Geometry::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Geometry::update_surface_area), "update_surface_area")
        
        
        .def("update", static_cast<void (PhysiCell::Geometry::*) (PhysiCell::Cell*, PhysiCell::Phenotype&, double)> (&PhysiCell::Geometry::update), "update")
        ;
        
    //Mechanics
    py::class_<PhysiCell::Mechanics>(pcore, "Mechanics")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("cell_cell_adhesion_strength", &PhysiCell::Mechanics::cell_cell_adhesion_strength)
        .def_readwrite("cell_BM_adhesion_strength", &PhysiCell::Mechanics::cell_BM_adhesion_strength)
        .def_readwrite("cell_cell_repulsion_strength", &PhysiCell::Mechanics::cell_cell_repulsion_strength)
        .def_readwrite("cell_BM_repulsion_strength", &PhysiCell::Mechanics::cell_BM_repulsion_strength)
        
        .def_readwrite("relative_maximum_adhesion_distance", &PhysiCell::Mechanics::relative_maximum_adhesion_distance)
        
        .def_readwrite("relative_maximum_attachment_distance", &PhysiCell::Mechanics::relative_maximum_attachment_distance)
        .def_readwrite("relative_detachment_distance", &PhysiCell::Mechanics::relative_detachment_distance)
        
        .def_readwrite("maximum_number_of_attachments", &PhysiCell::Mechanics::maximum_number_of_attachments)
        .def_readwrite("attachment_elastic_constant", &PhysiCell::Mechanics::attachment_elastic_constant)
        .def_readwrite("maximum_attachment_rate", &PhysiCell::Mechanics::maximum_attachment_rate)
        
        .def("set_relative_maximum_adhesion_distance", static_cast<void (PhysiCell::Mechanics::*) (double)> (&PhysiCell::Mechanics::set_relative_maximum_adhesion_distance), "set_relative_maximum_adhesion_distance")
        
        .def("set_relative_equilibrium_distance", static_cast<void (PhysiCell::Mechanics::*) (double)> (&PhysiCell::Mechanics::set_relative_equilibrium_distance), "set_relative_equilibrium_distance")
        
        
        .def("set_absolute_equilibrium_distance", static_cast<void (PhysiCell::Mechanics::*) (PhysiCell::Phenotype& , double)> (&PhysiCell::Mechanics::set_absolute_equilibrium_distance), "set_absolute_equilibrium_distance")
        
        ;
        
        
    //Motility
    py::class_<PhysiCell::Motility>(pcore, "Motility")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("is_motile", &PhysiCell::Motility::is_motile)
        .def_readwrite("persistence_time", &PhysiCell::Motility::persistence_time)
        .def_readwrite("migration_speed", &PhysiCell::Motility::migration_speed)
        .def_readwrite("migration_bias_direction", &PhysiCell::Motility::migration_bias_direction)
        .def_readwrite("migration_bias", &PhysiCell::Motility::migration_bias)
        .def_readwrite("restrict_to_2D", &PhysiCell::Motility::restrict_to_2D)
        .def_readwrite("motility_vector", &PhysiCell::Motility::motility_vector)
        .def_readwrite("chemotaxis_index", &PhysiCell::Motility::chemotaxis_index)
        .def_readwrite("chemotaxis_direction", &PhysiCell::Motility::chemotaxis_direction)
        ;
        
    //Secretion
    py::class_<PhysiCell::Secretion>(pcore, "Secretion")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("pMicroenvironment", &PhysiCell::Secretion::pMicroenvironment)
        
        .def_readwrite("secretion_rates", &PhysiCell::Secretion::secretion_rates)
        .def_readwrite("uptake_rates", &PhysiCell::Secretion::uptake_rates)
        .def_readwrite("saturation_densities", &PhysiCell::Secretion::saturation_densities)
        .def_readwrite("net_export_rates", &PhysiCell::Secretion::net_export_rates)
        
        //operators
        .def("sync_to_current_microenvironment", static_cast<void (PhysiCell::Secretion::*) ( void)> (&PhysiCell::Secretion::sync_to_current_microenvironment), "sync_to_current_microenvironment")
        
        .def("sync_to_microenvironment", static_cast<void (PhysiCell::Secretion::*) ( BioFVM::Microenvironment*)> (&PhysiCell::Secretion::sync_to_microenvironment), "sync_to_microenvironment")
        
        .def("advance", static_cast<void (PhysiCell::Secretion::*) ( BioFVM::Basic_Agent*, PhysiCell::Phenotype&, double)> (&PhysiCell::Secretion::advance), "advance")
        
        .def("set_all_secretion_to_zero", static_cast<void (PhysiCell::Secretion::*) (void)> (&PhysiCell::Secretion::set_all_secretion_to_zero), "set_all_secretion_to_zero")
        
        .def("set_all_uptake_to_zero", static_cast<void (PhysiCell::Secretion::*) (void)> (&PhysiCell::Secretion::set_all_uptake_to_zero), "set_all_uptake_to_zero")
        
        .def("scale_all_secretion_by_factor", static_cast<void (PhysiCell::Secretion::*) (double)> (&PhysiCell::Secretion::scale_all_secretion_by_factor), "scale_all_secretion_by_factor")
        
        .def("scale_all_uptake_by_factor", static_cast<void (PhysiCell::Secretion::*) (double)> (&PhysiCell::Secretion::scale_all_uptake_by_factor), "scale_all_uptake_by_factor")
        
        ;
        
    //Cell_Functions
    py::class_<PhysiCell::Cell_Functions>(pcore, "Cell_Functions")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("cycle_model", &PhysiCell::Cell_Functions::cycle_model)
        
        
        //TODO: Function pointers
        //void (*volume_update_function)( Cell* pCell, Phenotype& phenotype , double dt );
        //.def("volume_update_function", static_cast<void (PhysiCell::CellFunctions::*) (PhysiCell::Cell*, PhysiCell::Phenotype&,double)> (&PhysiCell::Cell_Functions::volume_update_function))
        //.def("assign_volume_update_function", [](Physicell::Cell_Functions &cellfunction){cellfunction.volume_update_function})
        ;
    //Molecular
    py::class_<PhysiCell::Molecular>(pcore, "Molecular")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("internalized_total_substrates", &PhysiCell::Molecular::internalized_total_substrates)
        
        .def_readwrite("fraction_released_at_death", &PhysiCell::Molecular::fraction_released_at_death)
        
        .def_readwrite("fraction_transferred_when_ingested", &PhysiCell::Molecular::fraction_transferred_when_ingested)
        
        //operators
        .def("sync_to_current_microenvironment", static_cast<void (PhysiCell::Molecular::*) ( void)> (&PhysiCell::Molecular::sync_to_current_microenvironment), "sync_to_current_microenvironment")
        
        .def("sync_to_microenvironment", static_cast<void (PhysiCell::Molecular::*) ( BioFVM::Microenvironment*)> (&PhysiCell::Molecular::sync_to_microenvironment), "sync_to_microenvironment")
        
        .def("sync_to_cell", static_cast<void (PhysiCell::Molecular::*) ( BioFVM::Basic_Agent*)> (&PhysiCell::Molecular::sync_to_cell), "sync_to_cell")
        
        ;
    //Phenotype
    py::class_<PhysiCell::Phenotype>(pcore, "Phenotype")
        //contructors
        .def(py::init<>())
        
        //attributes
        .def_readwrite("flagged_for_division", &PhysiCell::Phenotype::flagged_for_division)
        .def_readwrite("flagged_for_removal", &PhysiCell::Phenotype::flagged_for_removal)
        
        .def_readwrite("cycle", &PhysiCell::Phenotype::cycle)
        .def_readwrite("death", &PhysiCell::Phenotype::death)
        .def_readwrite("volume", &PhysiCell::Phenotype::volume)
        .def_readwrite("geometry", &PhysiCell::Phenotype::geometry)
        .def_readwrite("mechanics", &PhysiCell::Phenotype::mechanics)
        .def_readwrite("motility", &PhysiCell::Phenotype::motility)
        .def_readwrite("secretion", &PhysiCell::Phenotype::secretion)
        
        .def_readwrite("molecular", &PhysiCell::Phenotype::molecular)
        
        //operators
        .def("sync_to_functions", static_cast<void (PhysiCell::Phenotype::*) ( PhysiCell::Cell_Functions&)> (&PhysiCell::Phenotype::sync_to_functions), "sync_to_functions")
        
        .def("sync_to_microenvironment", static_cast<void (PhysiCell::Phenotype::*) ( BioFVM::Microenvironment*)> (&PhysiCell::Phenotype::sync_to_microenvironment), "sync_to_microenvironment")
        
        //.def("sync_to_default_functions", static_cast<void (PhysiCell::Phenotype::*) ( void)> (&PhysiCell::Phenotype::sync_to_default_functions), "sync_to_default_functions")
        ;
        
    //Objects and Functions Designed for use with Python
        
    //Cell Container_py (replacing the legacy cell container to better maintain the update_all_cells
    py::class_<PhysiCellCore_py::Cell_Container_py>(pcore, "Cell_Container")
    //contructors
        .def(py::init<>())
        .def(py::init<BioFVM::Microenvironment& , double>(), py::arg("microenvironment"), py::arg("voxel_size"))
        
        //attributes
        .def_readwrite("all_cells", &PhysiCellCore_py::Cell_Container_py::all_cells)
        
        //operators
        .def("initialize", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (double, double, double, double, double, double, double)> (&PhysiCellCore_py::Cell_Container_py::initialize), "initialize", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("voxel_size"))
        .def("initialize", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (double, double, double, double, double, double, double, double, double)> (&PhysiCellCore_py::Cell_Container_py::initialize), "initialize", py::arg("x_start"), py::arg("x_end"), py::arg("y_start"), py::arg("y_end"), py::arg("z_start"), py::arg("z_end"), py::arg("dx"), py::arg("dy"), py::arg("dz"))
        
        .def("create_cell_container_for_microenvironment", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (BioFVM::Microenvironment&, double)> (&PhysiCellCore_py::Cell_Container_py::create_cell_container_for_microenvironment), "create_cell_container_for_microenvironment", py::arg("microenvironment"), py::arg("voxel_size"))
        
        .def("update_all_cells", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (double, double, double)> (&PhysiCellCore_py::Cell_Container_py::update_all_cells), "update_all_cells", py::arg("t"), py::arg("phenotype_dt"), py::arg("mechanics_dt"))
        .def("update_all_cells", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (double, double, double, double)> (&PhysiCellCore_py::Cell_Container_py::update_all_cells), "update_all_cells", py::arg("t"), py::arg("phenotype_dt"), py::arg("mechanics_dt"), py::arg("diffusion_dt"))
        
        .def("create_cell", static_cast<PhysiCell::Cell* (PhysiCellCore_py::Cell_Container_py::*) (void)> (&PhysiCellCore_py::Cell_Container_py::create_cell), "create_cell")
        .def("create_cell", static_cast<PhysiCell::Cell* (PhysiCellCore_py::Cell_Container_py::*) (PhysiCell::Cell_Definition&)> (&PhysiCellCore_py::Cell_Container_py::create_cell), "create_cell")
        
        .def("delete_cell", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (int)> (&PhysiCellCore_py::Cell_Container_py::delete_cell), "delete_cell", py::arg("cell_index"))
        .def("delete_cell", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (PhysiCell::Cell*)> (&PhysiCellCore_py::Cell_Container_py::delete_cell), "delete_cell", py::arg("cell"))
        
        .def("delete_cell_original", static_cast<void (PhysiCellCore_py::Cell_Container_py::*) (int)> (&PhysiCellCore_py::Cell_Container_py::delete_cell_original), "delete_cell_original", py::arg("cell_index"))
        ;
    pcore.def("save_MultiCellDS", &PhysiCellCore_py::save_PhysiCell_to_MultiCellDS_xml_pugi_py, "Save Physicell multicell datastructure", py::arg("filename"), py::arg("microenvironment"), py::arg("cell_container"), py::arg("current_simulation_time"));
    
    //pcore.def("create_cell_in_cellcontainer", &PhysiCellCore_py::create_cell_py , "add cells to the cell container", py::arg("cell container"));

};
