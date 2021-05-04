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

#include "biofvm_py.h"


namespace BioFVM_py
{

Microenvironment_py::Microenvironment_py()
{
    name = "unnamed"; 
    spatial_units = "none"; 
    time_units = "none";

    bulk_source_sink_solver_setup_done = false; 
    thomas_setup_done = false; 
    diffusion_solver_setup_done = false; 

    //diffusion_decay_solver = empty_diffusion_solver_py;
    diffusion_decay_solver = BioFVM::diffusion_decay_solver__constant_coefficients_LOD_3D; 

    mesh.resize(1,1,1); 

    one.resize( 1 , 1.0 ); 
    zero.resize( 1 , 0.0 );

    temporary_density_vectors1.resize( mesh.voxels.size() , zero ); 
    temporary_density_vectors2.resize( mesh.voxels.size() , zero ); 
    p_density_vectors = &temporary_density_vectors1;

    gradient_vectors.resize( mesh.voxels.size() ); 
    for( unsigned int k=0 ; k < mesh.voxels.size() ; k++ )
    {
        gradient_vectors[k].resize( 1 ); 
        (gradient_vectors[k])[0].resize( 3, 0.0 );
    }
    gradient_vector_computed.resize( mesh.voxels.size() , false ); 

    bulk_supply_rate_function = zero_function; 
    bulk_supply_target_densities_function = zero_function; 
    bulk_uptake_rate_function = zero_function; 

    density_names.assign( 1 , "unnamed" ); 
    density_units.assign( 1 , "none" ); 

    diffusion_coefficients.assign( number_of_densities() , 0.0 ); 
    decay_rates.assign( number_of_densities() , 0.0 ); 

    one_half = one; 
    one_half *= 0.5; 

    one_third = one; 
    one_third /= 3.0;

    dirichlet_value_vectors.assign( mesh.voxels.size(), one ); 
    dirichlet_activation_vector.assign( 1 , true ); 

    dirichlet_activation_vectors.assign( 1 , dirichlet_activation_vector ); 

    microenvironment_options = default_microenvironment_options;
    
    microenvironment_options.Dirichlet_all.assign( 1 , true ); 
    microenvironment_options.Dirichlet_xmin.assign( 1 , false ); 
    microenvironment_options.Dirichlet_xmax.assign( 1 , false ); 
    microenvironment_options.Dirichlet_ymin.assign( 1 , false ); 
    microenvironment_options.Dirichlet_ymax.assign( 1 , false ); 
    microenvironment_options.Dirichlet_zmin.assign( 1 , false ); 
    microenvironment_options.Dirichlet_zmax.assign( 1 , false ); 

    microenvironment_options.Dirichlet_xmin_values.assign( 1 , 1.0 ); 
    microenvironment_options.Dirichlet_xmax_values.assign( 1 , 1.0 ); 
    microenvironment_options.Dirichlet_ymin_values.assign( 1 , 1.0 ); 
    microenvironment_options.Dirichlet_ymax_values.assign( 1 , 1.0 ); 
    microenvironment_options.Dirichlet_zmin_values.assign( 1 , 1.0 ); 
    microenvironment_options.Dirichlet_zmax_values.assign( 1 , 1.0 ); 


    return;
};

Microenvironment_py::Microenvironment_py(BioFVM::Microenvironment_Options microenvironment_options_in)
{
    Microenvironment_py();
    microenvironment_options = microenvironment_options_in;
    return;
};
 

Microenvironment_py::Microenvironment_py(std::string filename)
{

    Microenvironment_py(root_node);
    setup_microenvironment_from_file(filename);
    return;
};
Microenvironment_py::Microenvironment_py(pugi::xml_node root_node )
{
    Microenvironment_py();
    setup_microenvironment_from_XML_node( root_node );
    return;
};


bool Microenvironment_py::setup_microenvironment_from_XML_node( pugi::xml_node root_node )
{
    pugi::xml_node node; 

    // First, look for the correct XML node. 
    // If it isn't there, return false. 

    node = PhysiCell::xml_find_node( root_node , "microenvironment_setup" );
    if( !node )
    { return false; }

    // now that we're using the XML to specify the microenvironment, don't 
    // use old defaults 

    // Don't let BioFVM use oxygen as the default 

    microenvironment_options.use_oxygen_as_first_field = false; 

    std::vector<double> initial_condition_vector = {}; 
    std::vector<double> Dirichlet_condition_vector = {}; 
    std::vector<bool> Dirichlet_activation_vector = {}; 

    std::vector<bool> Dirichlet_all = {}; 
    std::vector<bool> Dirichlet_xmin = {}; 
    std::vector<bool> Dirichlet_xmax = {}; 
    std::vector<bool> Dirichlet_ymin = {}; 
    std::vector<bool> Dirichlet_ymax = {}; 
    std::vector<bool> Dirichlet_zmin = {}; 
    std::vector<bool> Dirichlet_zmax = {}; 

    std::vector<double> Dirichlet_xmin_values = {}; 
    std::vector<double> Dirichlet_xmax_values = {}; 
    std::vector<double> Dirichlet_ymin_values = {}; 
    std::vector<double> Dirichlet_ymax_values = {}; 
    std::vector<double> Dirichlet_zmin_values = {}; 
    std::vector<double> Dirichlet_zmax_values = {}; 
    std::vector<double> Dirichlet_interior_values = {}; 


    // next, add all the substrates to the microenvironment
    // build the initial conditions and Dirichlet conditions as we go 

    // find the first substrate 
    pugi::xml_node node1 = node.child( "variable" ); // xml_find_node( node , "variable" ); 
    node = node1; 
    int i = 0; 

    bool activated_Dirichlet_boundary_detected = false; 

    while( node )
    {
        // get the name and units 
        std::string name = node.attribute( "name" ).value(); 
        std::string units = node.attribute( "units" ).value(); 
        
        // add the substrate 
        if( i == 0 )
        { set_density( 0, name, units ); }
        else
        { add_density( name, units ); }
        
        // get the diffusion and decay parameters 
        node1 = PhysiCell::xml_find_node( node, "physical_parameter_set" ); 
        
        diffusion_coefficients[i] = 
            PhysiCell::xml_get_double_value( node1, "diffusion_coefficient" ); 
        decay_rates[i] = 
            PhysiCell::xml_get_double_value( node1, "decay_rate" ); 
            
        // now, get the initial value  
        node1 = PhysiCell::xml_find_node( node, "initial_condition" ); 
        initial_condition_vector.push_back( PhysiCell::xml_get_my_double_value(node1) );
        
        // now, get the Dirichlet value
        node1 = PhysiCell::xml_find_node( node, "Dirichlet_boundary_condition" ); 
        Dirichlet_condition_vector.push_back( PhysiCell::xml_get_my_double_value(node1) );

        // now, decide whether or not to enable it 
        Dirichlet_activation_vector.push_back( node1.attribute("enabled").as_bool() );

        Dirichlet_all.push_back( Dirichlet_activation_vector[i] ); 
        if( Dirichlet_activation_vector[i] )
        { activated_Dirichlet_boundary_detected = true; }
        
        // default interior activation will mirror the boundary 
        
        Dirichlet_xmin.push_back( Dirichlet_activation_vector[i] ); 
        Dirichlet_xmax.push_back( Dirichlet_activation_vector[i] ); 
        Dirichlet_ymin.push_back( Dirichlet_activation_vector[i] ); 
        Dirichlet_ymax.push_back( Dirichlet_activation_vector[i] ); 
        Dirichlet_zmin.push_back( Dirichlet_activation_vector[i] ); 
        Dirichlet_zmax.push_back( Dirichlet_activation_vector[i] ); 
        
        Dirichlet_xmin_values.push_back( Dirichlet_condition_vector[i] ); 
        Dirichlet_xmax_values.push_back( Dirichlet_condition_vector[i] ); 
        Dirichlet_ymin_values.push_back( Dirichlet_condition_vector[i] ); 
        Dirichlet_ymax_values.push_back( Dirichlet_condition_vector[i] ); 
        Dirichlet_zmin_values.push_back( Dirichlet_condition_vector[i] ); 
        Dirichlet_zmax_values.push_back( Dirichlet_condition_vector[i] ); 
        
        // now figure out finer-grained controls 
        
        node1 = node.child( "Dirichlet_options" );
        if( node1 )
        {
            // xmin, xmax, ymin, ymax, zmin, zmax, interior 
            pugi::xml_node node2 = node1.child("boundary_value"); 
            
            while( node2 )
            {
                // which boundary? 
                std::string boundary_ID = node2.attribute("ID").value(); 
                
                // xmin 
                if( std::strstr( boundary_ID.c_str() , "xmin" ) )
                {
                    // on or off 
                    Dirichlet_xmin[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }
                    
                    // which value 
                    { Dirichlet_xmin_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }
                
                // xmax 
                if( std::strstr( boundary_ID.c_str() , "xmax" ) )
                {
                    // on or off 
                    Dirichlet_xmax[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }
                
                    // which value 
                    { Dirichlet_xmax_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }
                
                // ymin 
                if( std::strstr( boundary_ID.c_str() , "ymin" ) )
                {
                    // on or off 
                    Dirichlet_ymin[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }
                
                    // which value 
                    { Dirichlet_ymin_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }
                
                // ymax 
                if( std::strstr( boundary_ID.c_str() , "ymax" ) )
                {
                    // on or off 
                    Dirichlet_ymax[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }					
                    
                    // which value 
                    { Dirichlet_ymax_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }				
                                
                // zmin 
                if( std::strstr( boundary_ID.c_str() , "zmin" ) )
                {
                    // on or off 
                    Dirichlet_zmin[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }
                
                    // which value 
                    { Dirichlet_zmin_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }
                
                // zmax 
                if( std::strstr( boundary_ID.c_str() , "zmax" ) )
                {
                    // on or off 
                    Dirichlet_zmax[i] = node2.attribute("enabled").as_bool();
                    // if there is at least one off bondary here, "all" is false for this substrate 
                    if( node2.attribute("enabled").as_bool() == false )
                    { Dirichlet_all[i] = false; }
                
                    // which value 
                    { Dirichlet_zmax_values[i] = PhysiCell::xml_get_my_double_value( node2 ); }
                }
                
                node2 = node2.next_sibling("boundary_value"); 
            }
        }
        
        
        // move on to the next variable (if any!)
        node = node.next_sibling( "variable" ); 
        i++; 
    }

    // now that all the variables and boundary / initial conditions are defined, 
    // make sure that BioFVM knows about them 

    microenvironment_options.Dirichlet_condition_vector = Dirichlet_condition_vector;  
    microenvironment_options.Dirichlet_activation_vector = Dirichlet_activation_vector;
    microenvironment_options.initial_condition_vector = initial_condition_vector; 

    microenvironment_options.Dirichlet_all = Dirichlet_all; 

    microenvironment_options.Dirichlet_xmin = Dirichlet_xmin; 
    microenvironment_options.Dirichlet_xmax = Dirichlet_xmax; 
    microenvironment_options.Dirichlet_ymin = Dirichlet_ymin; 
    microenvironment_options.Dirichlet_ymax = Dirichlet_ymax; 
    microenvironment_options.Dirichlet_zmin = Dirichlet_zmin; 
    microenvironment_options.Dirichlet_zmax = Dirichlet_zmax; 

    microenvironment_options.Dirichlet_xmin_values = Dirichlet_xmin_values; 
    microenvironment_options.Dirichlet_xmax_values = Dirichlet_xmax_values; 
    microenvironment_options.Dirichlet_ymin_values = Dirichlet_ymin_values; 
    microenvironment_options.Dirichlet_ymax_values = Dirichlet_ymax_values; 
    microenvironment_options.Dirichlet_zmin_values = Dirichlet_zmin_values; 
    microenvironment_options.Dirichlet_zmax_values = Dirichlet_zmax_values; 

    // because outer boundary Dirichlet conditions are defined in the XML, 
    // make sure we don't accidentally disable them 

    microenvironment_options.outer_Dirichlet_conditions = false;

    // if *any* of the substrates have outer Dirichlet conditions enables, 
    // then set teh outer_Dirichlet_conditions = true; 

    if( activated_Dirichlet_boundary_detected ) 
    {
        microenvironment_options.outer_Dirichlet_conditions = true;
    }

    std::cout << activated_Dirichlet_boundary_detected << std::endl; 
    std::cout << "dc? " << microenvironment_options.outer_Dirichlet_conditions << std::endl; 

    // now, get the options 
    node = PhysiCell::xml_find_node( root_node , "microenvironment_setup" );
    node = PhysiCell::xml_find_node( node , "options" ); 

    // calculate gradients? 
    microenvironment_options.calculate_gradients = PhysiCell::xml_get_bool_value( node, "calculate_gradients" ); 

    // track internalized substrates in each agent? 
    microenvironment_options.track_internalized_substrates_in_each_agent 
        = PhysiCell::xml_get_bool_value( node, "track_internalized_substrates_in_each_agent" ); 

    // not yet supported : read initial conditions 
    /*
    // read in initial conditions from an external file 
            <!-- not yet supported --> 
            <initial_condition type="matlab" enabled="false">
                <filename>./config/initial.mat</filename>
            </initial_condition>
    */

    // not yet supported : read Dirichlet nodes (including boundary)
    /*
    // Read in Dirichlet nodes from an external file.
    // Note that if they are defined this way, then 
    // set 	default_microenvironment_options.outer_Dirichlet_conditions = false;
    // so that the microenvironment initialization in BioFVM does not 
    // also add Dirichlet nodes at the outer boundary

            <!-- not yet supported --> 
            <dirichlet_nodes type="matlab" enabled="false">
                <filename>./config/dirichlet.mat</filename>
            </dirichlet_nodes>
    */	

    
    node = PhysiCell::xml_find_node( root_node , "domain" );

    double xmin = PhysiCell::xml_get_double_value( node , "x_min" );
    double xmax = PhysiCell::xml_get_double_value( node , "x_max" );
    double ymin = PhysiCell::xml_get_double_value( node , "y_min" );
    double ymax = PhysiCell::xml_get_double_value( node , "y_max" );
    double zmin = PhysiCell::xml_get_double_value( node , "z_min" );
    double zmax = PhysiCell::xml_get_double_value( node , "z_max" );
    double dx = PhysiCell::xml_get_double_value( node, "dx" ); 
    double dy = PhysiCell::xml_get_double_value( node, "dy" ); 
    double dz = PhysiCell::xml_get_double_value( node, "dz" ); 

    microenvironment_options.simulate_2D = PhysiCell::xml_get_bool_value( node, "use_2D" ); 

    if( microenvironment_options.simulate_2D == true )
    {
        zmin = -0.5 * dz; 
        zmax = 0.5 * dz; 
    }			
    microenvironment_options.X_range = {xmin, xmax}; 
    microenvironment_options.Y_range = {ymin, ymax}; 
    microenvironment_options.Z_range = {zmin, zmax}; 

    microenvironment_options.dx = dx; 
    microenvironment_options.dy = dy; 
    microenvironment_options.dz = dz;


    return true;  
};

void Microenvironment_py::setup_microenvironment_from_file( std::string filename )
{
    pugi::xml_document physicell_config_doc;
    pugi::xml_node physicell_config_root; 
    
    std::cout << "Using config file for microenvironment:" << filename << " ... " << std::endl ; 
    pugi::xml_parse_result result = physicell_config_doc.load_file( filename.c_str()  );

    if( result.status != pugi::xml_parse_status::status_ok )
    {
        std::cout << "Error loading " << filename << "!" << std::endl; 
        return;
    }

    physicell_config_root = physicell_config_doc.child("PhysiCell_settings");
    setup_microenvironment_from_XML_node(physicell_config_root);
    return;
    return;
}

/*Add Density Functions*/

void Microenvironment_py::add_density_helper( void )
{
    // update sources and such 
	for( unsigned int i=0; i < temporary_density_vectors1.size() ; i++ )
	{
		temporary_density_vectors1[i].push_back( 0.0 ); 
		temporary_density_vectors2[i].push_back( 0.0 ); 
	}

	// resize the gradient data structures 
	for( unsigned int k=0 ; k < mesh.voxels.size() ; k++ )
	{
		gradient_vectors[k].resize( number_of_densities() ); 
		for( unsigned int i=0 ; i < number_of_densities() ; i++ )
		{
			(gradient_vectors[k])[i].resize( 3, 0.0 );
		}
	}

	gradient_vector_computed.resize( mesh.voxels.size() , false ); 	
	
	one_half = one; 
	one_half *= 0.5; 
	
	one_third = one; 
	one_third /= 3.0; 
	
	dirichlet_value_vectors.assign( mesh.voxels.size(), one ); 
	dirichlet_activation_vector.push_back( true ); 
	dirichlet_activation_vectors.assign( mesh.voxels.size(), dirichlet_activation_vector ); 
	
	// Fixes in PhysiCell preview November 2017
	microenvironment_options.Dirichlet_condition_vector.push_back( 1.0 ); //  = one; 
	microenvironment_options.Dirichlet_activation_vector.push_back( true ); // assign( number_of_densities(), true ); 
	
	microenvironment_options.initial_condition_vector.push_back( 1.0 ); 

	microenvironment_options.Dirichlet_all.push_back( true ); 
//	default_microenvironment_options.Dirichlet_interior.push_back( true );
	microenvironment_options.Dirichlet_xmin.push_back( false ); 
	microenvironment_options.Dirichlet_xmax.push_back( false ); 
	microenvironment_options.Dirichlet_ymin.push_back( false ); 
	microenvironment_options.Dirichlet_ymax.push_back( false ); 
	microenvironment_options.Dirichlet_zmin.push_back( false ); 
	microenvironment_options.Dirichlet_zmax.push_back( false ); 
	
	microenvironment_options.Dirichlet_xmin_values.push_back( 1.0 ); 
	microenvironment_options.Dirichlet_xmax_values.push_back( 1.0 ); 
	microenvironment_options.Dirichlet_ymin_values.push_back( 1.0 ); 
	microenvironment_options.Dirichlet_ymax_values.push_back( 1.0 ); 
	microenvironment_options.Dirichlet_zmin_values.push_back( 1.0 ); 
	microenvironment_options.Dirichlet_zmax_values.push_back( 1.0 ); 
	
    return;
    
}

void Microenvironment_py::resize_densities( int new_size )
{
	zero.assign( new_size, 0.0 ); 
	one.assign( new_size , 1.0 );

	temporary_density_vectors1.assign( mesh.voxels.size() , zero );
	temporary_density_vectors2.assign( mesh.voxels.size() , zero );

	for( unsigned int k=0 ; k < mesh.voxels.size() ; k++ )
	{
		gradient_vectors[k].resize( number_of_densities() ); 
		for( unsigned int i=0 ; i < number_of_densities() ; i++ )
		{
			(gradient_vectors[k])[i].resize( 3, 0.0 );
		}
	}
	gradient_vector_computed.resize( mesh.voxels.size() , false ); 	
	
	diffusion_coefficients.assign( new_size , 0.0 ); 
	decay_rates.assign( new_size , 0.0 ); 

	density_names.assign( new_size, "unnamed" ); 
	density_units.assign( new_size , "none" ); 

	one_half = one; 
	one_half *= 0.5; 
	
	one_third = one; 
	one_third /= 3.0; 
	
	dirichlet_value_vectors.assign( mesh.voxels.size(), one ); 
	dirichlet_activation_vector.assign( new_size, true ); 

	dirichlet_activation_vectors.assign( mesh.voxels.size(), dirichlet_activation_vector ); 

	microenvironment_options.Dirichlet_condition_vector.assign( new_size , 1.0 );  
	microenvironment_options.Dirichlet_activation_vector.assign( new_size, true ); 
	
	microenvironment_options.initial_condition_vector.assign( new_size , 1.0 ); 
	
	microenvironment_options.Dirichlet_all.assign( new_size , true ); 
//	default_microenvironment_options.Dirichlet_interior.assign( new_size, true );
	microenvironment_options.Dirichlet_xmin.assign( new_size , false ); 
	microenvironment_options.Dirichlet_xmax.assign( new_size , false ); 
	microenvironment_options.Dirichlet_ymin.assign( new_size , false ); 
	microenvironment_options.Dirichlet_ymax.assign( new_size , false ); 
	microenvironment_options.Dirichlet_zmin.assign( new_size , false ); 
	microenvironment_options.Dirichlet_zmax.assign( new_size , false ); 
	
	microenvironment_options.Dirichlet_xmin_values.assign( new_size , 1.0 ); 
	microenvironment_options.Dirichlet_xmax_values.assign( new_size , 1.0 ); 
	microenvironment_options.Dirichlet_ymin_values.assign( new_size , 1.0 ); 
	microenvironment_options.Dirichlet_ymax_values.assign( new_size , 1.0 ); 
	microenvironment_options.Dirichlet_zmin_values.assign( new_size , 1.0 ); 
	microenvironment_options.Dirichlet_zmax_values.assign( new_size , 1.0 ); 

	return; 
}


void Microenvironment_py::add_density( void )
{
	// fix in PhysiCell preview November 2017 
	// default_microenvironment_options.use_oxygen_as_first_field = false; 
	
	// update 1, 0 
	zero.push_back( 0.0 ); 
	one.push_back( 1.0 );
	
	// update units
	density_names.push_back( "unnamed" ); 
	density_units.push_back( "none" ); 

	// update coefficients 
	diffusion_coefficients.push_back( 0.0 ); 
	decay_rates.push_back( 0.0 ); 
	
    add_density_helper();
	
	return; 
}

void Microenvironment_py::add_density( std::string name , std::string units )
{
	// fix in PhysiCell preview November 2017 
	// default_microenvironment_options.use_oxygen_as_first_field = false; 
	
	// update 1, 0 
	zero.push_back( 0.0 ); 
	one.push_back( 1.0 );

	// update units
	density_names.push_back( name ); 
	density_units.push_back( units ); 

	// update coefficients 
	diffusion_coefficients.push_back( 0.0 ); 
	decay_rates.push_back( 0.0 ); 
	
	add_density_helper(); 
	
	return; 
}

void Microenvironment_py::add_density( std::string name , std::string units, double diffusion_constant, double decay_rate )
{
	// fix in PhysiCell preview November 2017 
	// default_microenvironment_options.use_oxygen_as_first_field = false; 
	
	// update 1, 0 
	zero.push_back( 0.0 ); 
	one.push_back( 1.0 );
	
	// update units
	density_names.push_back( name ); 
	density_units.push_back( units ); 

	// update coefficients 
	diffusion_coefficients.push_back( diffusion_constant ); 
	decay_rates.push_back( decay_rate );
 
    add_density_helper(); 
    
    return;
}

void Microenvironment_py::set_density( int index , std::string name , std::string units )
{
	// fix in PhysiCell preview November 2017 
	if( index == 0 )
	{ microenvironment_options.use_oxygen_as_first_field = false; }
	
	density_names[index] = name; 
	density_units[index] = units; 
	return; 
}

void Microenvironment_py::set_density( int index , std::string name , std::string units , double diffusion_constant , double decay_rate )
{
	// fix in PhysiCell preview November 2017 
	if( index == 0 )
	{ microenvironment_options.use_oxygen_as_first_field = false; }
	
	density_names[index] = name; 
	density_units[index] = units; 
	
	diffusion_coefficients[index] = diffusion_constant; 
	decay_rates[index] = decay_rate;	
	return; 
}


void Microenvironment_py::initialize_microenvironment( void )
{
	// create and name a microenvironment; 
	name = microenvironment_options.name;
	// register the diffusion solver 
	if( microenvironment_options.simulate_2D == true )
	{
		diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_2D; 
	}
	else
	{
		diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	}
	
	// set the default substrate to oxygen (with typical units of mmHg)
	if( microenvironment_options.use_oxygen_as_first_field == true )
	{
		set_density(0, "oxygen" , "mmHg" );
		diffusion_coefficients[0] = 1e5; 
		decay_rates[0] = 0.1; 
	}
	
	// resize the microenvironment  
	if( microenvironment_options.simulate_2D == true )
	{
		microenvironment_options.Z_range[0] = -microenvironment_options.dz/2.0; 
		microenvironment_options.Z_range[1] = microenvironment_options.dz/2.0;
	}
	resize_space( microenvironment_options.X_range[0], microenvironment_options.X_range[1] , 
		microenvironment_options.Y_range[0], microenvironment_options.Y_range[1], 
		microenvironment_options.Z_range[0], microenvironment_options.Z_range[1], 
		microenvironment_options.dx,microenvironment_options.dy,microenvironment_options.dz );
		
	// set units
	spatial_units = microenvironment_options.spatial_units;
	time_units = microenvironment_options.time_units;
	mesh.units = microenvironment_options.spatial_units;

	// set the initial densities to the values set in the initial_condition_vector
	
	// if the initial condition vector has not been set, use the Dirichlet condition vector 
	if( microenvironment_options.initial_condition_vector.size() != 
		number_of_densities() )
	{
		std::cout << "BioFVM Warning: Initial conditions not set. " << std::endl 
				  << "                Using Dirichlet condition vector to set initial substrate values!" << std::endl 
				  << "                In the future, set microenvironment_options.initial_condition_vector." 
				  << std::endl << std::endl;  
		microenvironment_options.initial_condition_vector = microenvironment_options.Dirichlet_condition_vector; 
	}

	// set the initial condition 
	for( unsigned int n=0; n < number_of_voxels() ; n++ )
	{ density_vector(n) = microenvironment_options.initial_condition_vector; }

	// now, figure out which sides have BCs (for at least one substrate): 

	bool xmin = false; 
	bool xmax = false; 
	bool ymin = false; 
	bool ymax = false; 
	bool zmin = false; 
	bool zmax = false; 
	
	if( microenvironment_options.outer_Dirichlet_conditions == true )
	{
		for( int n=0 ; n < number_of_densities() ; n++ )
		{
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_xmin[n] )
				{ xmin = true; }
			
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_xmax[n] )
				{ xmax = true; }
			
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_ymin[n] )
				{ ymin = true; }
			
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_ymax[n] )
				{ ymax = true; }
				
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_zmin[n] )
				{ zmin = true; }
			
			if( microenvironment_options.Dirichlet_all[n] || 
				microenvironment_options.Dirichlet_zmax[n] )
				{ zmax = true; }
		}
		
		// add the Dirichlet nodes in the right places 
		
	}
	std::cout << "which boundaries?" << std::endl; 
	std::cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << std::endl; 

	// add the Dirichlet nodes in the right places 
	// now, go in and set the values 
	if( microenvironment_options.outer_Dirichlet_conditions == true ) 
	{
		// set xmin if xmin = true or all = true 
		if( xmin == true )
		{
			for( unsigned int k=0 ; k < mesh.z_coordinates.size() ; k++ )
			{
				int I = 0; 
				// set Dirichlet conditions along the xmin outer edges 
				for( unsigned int j=0 ; j < mesh.y_coordinates.size() ; j++ )
				{
					// set the value 
					add_dirichlet_node( voxel_index(I,j,k) , microenvironment_options.Dirichlet_xmin_values );
					
					// set the activation 
					set_substrate_dirichlet_activation( voxel_index(I,j,k) , 
					microenvironment_options.Dirichlet_xmin ); 
					
				}
			}
		}			
		
		// set xmax if xmax = true or all = true 
		if( xmax == true )
		{
			for( unsigned int k=0 ; k < mesh.z_coordinates.size() ; k++ )
			{
				int I = mesh.x_coordinates.size()-1;; 
				// set Dirichlet conditions along the xmax outer edges 
				for( unsigned int j=0 ; j < mesh.y_coordinates.size() ; j++ )
				{
					// set the values 
					add_dirichlet_node( voxel_index(I,j,k) , microenvironment_options.Dirichlet_xmax_values );
					
					// set the activation 
					set_substrate_dirichlet_activation( voxel_index(I,j,k) , 
					microenvironment_options.Dirichlet_xmax ); 
				}
			}
		}			
		
		// set ymin if ymin = true or all = true 
		if( ymin == true )
		{
			for( unsigned int k=0 ; k < mesh.z_coordinates.size() ; k++ )
			{
				int J = 0; // mesh.x_coordinates.size()-1;; 
				// set Dirichlet conditions along the ymin outer edges 
				for( unsigned int i=0 ; i < mesh.x_coordinates.size() ; i++ )
				{
					// set the values 
					add_dirichlet_node( voxel_index(i,J,k) , microenvironment_options.Dirichlet_ymin_values );
					
					// set the activation 
					set_substrate_dirichlet_activation( voxel_index(i,J,k) , 
					microenvironment_options.Dirichlet_ymin ); 
				}
			}
		}	
		
		// set ymzx if ymax = true or all = true; 
		if( ymax == true )
		{
			for( unsigned int k=0 ; k < mesh.z_coordinates.size() ; k++ )
			{
				int J = mesh.y_coordinates.size()-1;; 
				// set Dirichlet conditions along the ymin outer edges 
				for( unsigned int i=0 ; i < mesh.x_coordinates.size() ; i++ )
				{
					// set the value 
					add_dirichlet_node( voxel_index(i,J,k) , microenvironment_options.Dirichlet_ymax_values );
					
					// set the activation 
					set_substrate_dirichlet_activation( voxel_index(i,J,k) , 
					microenvironment_options.Dirichlet_ymax ); 
				}
			}
		}	
		
		// if not 2D:
		if( microenvironment_options.simulate_2D == false )
		{
			// set zmin if zmin = true or all = true 
			if( zmin == true )
			{
				for( unsigned int j=0 ; j < mesh.y_coordinates.size() ; j++ )
				{
					int K = 0; // mesh.z_coordinates.size()-1;; 
					// set Dirichlet conditions along the ymin outer edges 
					for( unsigned int i=0 ; i < mesh.x_coordinates.size() ; i++ )
					{
						// set the value 
						add_dirichlet_node( voxel_index(i,j,K) , microenvironment_options.Dirichlet_zmin_values );
					
						// set the activation 
						set_substrate_dirichlet_activation( voxel_index(i,j,K) , 
						microenvironment_options.Dirichlet_zmin ); 
					}
				}
			}				
			
			// set zmax if zmax = true or all = true 
			if( zmax == true )
			{
				for( unsigned int j=0 ; j < mesh.y_coordinates.size() ; j++ )
				{
					int K = mesh.z_coordinates.size()-1;; 
					// set Dirichlet conditions along the ymin outer edges 
					for( unsigned int i=0 ; i < mesh.x_coordinates.size() ; i++ )
					{
						// set the value 
						add_dirichlet_node( voxel_index(i,j,K) , microenvironment_options.Dirichlet_zmax_values );
						
						// set the activation 
						set_substrate_dirichlet_activation( voxel_index(i,j,K) , 
						microenvironment_options.Dirichlet_zmax ); 						
					}
				}
			}				
		}
		
	}
	
	// set the Dirichlet condition activation vector to match the microenvironment options 
	for( int i=0 ; i < microenvironment_options.Dirichlet_activation_vector.size(); i++ )
	{
		set_substrate_dirichlet_activation( i , microenvironment_options.Dirichlet_activation_vector[i] ); 
	}
	
	display_information( std::cout );
	return;
}


void Microenvironment_py::display_information( std::ostream& os )
{
	os << std::endl << "Microenvironment summary: " << name << ": " << std::endl; 
	mesh.display_information( os ); 
	os << "Densities: (" << number_of_densities() << " total)" << std::endl; 
	for( unsigned int i = 0 ; i < density_names.size() ; i++ )
	{
		os << "   " << density_names[i] << ":" << std::endl
		<< "     units: " << density_units[i] << std::endl 
		<< "     diffusion coefficient: " << diffusion_coefficients[i]  
			<< " " << spatial_units << "^2 / " << time_units << std::endl
		<< "     decay rate: " << decay_rates[i] 
			<< " " << time_units << "^-1" << std::endl 
		<< "     diffusion length scale: " << sqrt( diffusion_coefficients[i] / ( 1e-12 + decay_rates[i] ) ) 
			<< " " << spatial_units << std::endl 
		<< "     initial condition: " << microenvironment_options.initial_condition_vector[i] 
			<< " " << density_units[i] << std::endl 
		<< "     boundary condition: " << microenvironment_options.Dirichlet_condition_vector[i] 
			<< " " << density_units[i] << " (enabled: "; 
		if( dirichlet_activation_vector[i] == true )
		{ os << "true"; }
		else
		{ os << "false"; }
		os << ")" << std::endl; 
	}
	os << std::endl; 
	
	return; 
}


void Microenvironment_py::display_info ( void )
{
 display_information(std::cout);   
    
}

}; //end namespace BioFVM


