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

#include "physicellcore_py.h"


namespace PhysiCellCore_py{
    
Cell_Container_py::Cell_Container_py()
{
    all_cells = (std::vector<PhysiCell::Cell*> *) &all_basic_agents;	
	boundary_condition_for_pushed_out_agents= PhysiCell::PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	std::vector<PhysiCell::Cell*> cells_ready_to_divide;
	std::vector<PhysiCell::Cell*> cells_ready_to_die;
	
	return; 
    
}

Cell_Container_py::Cell_Container_py(BioFVM::Microenvironment& m , double mechanics_voxel_size )
{
    Cell_Container_py();
    create_cell_container_for_microenvironment( m , mechanics_voxel_size );
	return; 
    
}

void Cell_Container_py::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size){
    initialize(x_start, x_end, y_start, y_end, z_start, z_end , voxel_size, voxel_size, voxel_size);
    return;
    
}

void Cell_Container_py::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz){
    all_cells = (std::vector<PhysiCell::Cell*> *) &all_basic_agents;	
	boundary_condition_for_pushed_out_agents= PhysiCell::PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	std::vector<PhysiCell::Cell*> cells_ready_to_divide;
	std::vector<PhysiCell::Cell*> cells_ready_to_die;

	underlying_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end , dx, dy, dz);
	agent_grid.resize(underlying_mesh.voxels.size());
	max_cell_interactive_distance_in_voxel.resize(underlying_mesh.voxels.size(), 0.0);
	agents_in_outer_voxels.resize(6);
	
	return; 
}

void Cell_Container_py::create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size ){
    initialize( m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
        m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
        m.mesh.bounding_box[2], m.mesh.bounding_box[5],  mechanics_voxel_size );
    m.agent_container = (BioFVM::Agent_Container*) this; //TODO FIX THIS
    return;
}

void Cell_Container_py::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ ){
    update_all_cells(t, phenotype_dt_, mechanics_dt_ , PhysiCell::diffusion_dt );
    return; 
}
void Cell_Container_py::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ ){
    	// secretions and uptakes. Syncing with BioFVM is automated. 

	// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "secretion" << std::endl; 
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		if( (*all_cells)[i]->is_out_of_domain == false )
		{
			(*all_cells)[i]->phenotype.secretion.advance( (*all_cells)[i], (*all_cells)[i]->phenotype , diffusion_dt_ );
		}
	}
	
	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	static double phenotype_dt_tolerance = 0.001 * phenotype_dt_; 
	static double mechanics_dt_tolerance = 0.001 * mechanics_dt_; 
	
	if( fabs(time_since_last_cycle-phenotype_dt_ ) < phenotype_dt_tolerance || !initialzed)
	{
		// Reset the max_radius in each voxel. It will be filled in set_total_volume
		// It might be better if we calculate it before mechanics each time 
		// std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);
		
		if(!initialzed)
		{
			time_since_last_cycle = phenotype_dt_;
		}
		
		// new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update, 
		// checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)
		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "bundled phenotype" << std::endl; 
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if( (*all_cells)[i]->is_out_of_domain == false )
			{
				(*all_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle ); 
			}
		}
		
		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "divide / die " << std::endl; 
		// process divides / removes 
		for( int i=0; i < cells_ready_to_divide.size(); i++ )
		{
			divide_cell(cells_ready_to_divide[i]); //TODO
		}
		for( int i=0; i < cells_ready_to_die.size(); i++ )
		{	
			kill_cell(cells_ready_to_die[i]);	//TODO
		}
		num_divisions_in_current_step+=  cells_ready_to_divide.size();
		num_deaths_in_current_step+=  cells_ready_to_die.size();
		
		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
		last_cell_cycle_time= t;
	}
		
	double time_since_last_mechanics= t- last_mechanics_time;
	
	// if( time_since_last_mechanics>= mechanics_dt || !initialzed)
	if( fabs(time_since_last_mechanics - mechanics_dt_) < mechanics_dt_tolerance || !initialzed)
	{
		if(!initialzed)
		{
			time_since_last_mechanics = mechanics_dt_;
		}
		
		// new February 2018 
		// if we need gradients, compute them
		if( default_microenvironment_options.calculate_gradients ) 
		{ microenvironment.compute_all_gradient_vectors();  }
		// end of new in Feb 2018 
		
		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "interactions" << std::endl; 
		// perform interactions -- new in June 2020 
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			PhysiCell::Cell* pC = (*all_cells)[i]; 
			if( pC->functions.contact_function && pC->is_out_of_domain == false )
			{ evaluate_interactions( pC,pC->phenotype,time_since_last_mechanics ); }
		}
		
		// perform custom computations 

		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "custom" << std::endl; 
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			PhysiCell::Cell* pC = (*all_cells)[i]; 
			if( pC->functions.custom_cell_rule && pC->is_out_of_domain == false )
			{ pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics ); }
		}
		
		// update velocities 
		
		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "velocity" << std::endl; 
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			PhysiCell::Cell* pC = (*all_cells)[i]; 
			if( pC->functions.update_velocity && pC->is_out_of_domain == false && pC->is_movable )
			{ pC->functions.update_velocity( pC,pC->phenotype,time_since_last_mechanics ); }
		}

		// update positions 
		
		// std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " " << "position" << std::endl; 
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			PhysiCell::Cell* pC = (*all_cells)[i]; 
			if( pC->is_out_of_domain == false && pC->is_movable)
			{ pC->update_position(time_since_last_mechanics); }
		}
		
		for( int i=0; i < (*all_cells).size(); i++ )
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
				(*all_cells)[i]->update_voxel_in_container();
		last_mechanics_time=t;
	}
	
	initialzed=true;
	return;
}

// making the create and save cell functions part of the 
PhysiCell::Cell* Cell_Container_py::create_cell( void ){
	std::cout << "1";
	PhysiCell::Cell* pNew = new PhysiCell::Cell;	
	std::cout << "3";	 
    (*all_cells).push_back( pNew ); 
		std::cout << "2";
    pNew->index=(*all_cells).size()-1;
	std::cout << "4";

    // All the phenotype and other data structures are already set 
    // by virtue of the default Cell constructor. 
std::cout << "5";
    pNew->set_total_volume( pNew->phenotype.volume.total ); 

    return pNew; 
}

PhysiCell::Cell* create_cell_py( Cell_Container_py& cellcontainer ){
	std::cout << "1";
	PhysiCell::Cell* pNew;
	pNew = &PhysiCell::Cell;	
	std::cout << "3";	 
    (*cellcontainer.all_cells).push_back( pNew ); 
		std::cout << "2";
    pNew->index=(*cellcontainer.all_cells).size()-1;
	std::cout << "4";

    // All the phenotype and other data structures are already set 
    // by virtue of the default Cell constructor. 
std::cout << "5";
    pNew->set_total_volume( pNew->phenotype.volume.total ); 

    return pNew; 
}



PhysiCell::Cell* Cell_Container_py::create_cell( PhysiCell::Cell_Definition& cd ){
    PhysiCell::Cell* pNew = create_cell(); 
	
	// use the cell defaults; 
	pNew->type = cd.type; 
	pNew->type_name = cd.name; 
	
	pNew->custom_data = cd.custom_data; 
	pNew->parameters = cd.parameters; 
	pNew->functions = cd.functions; 
	
	pNew->phenotype = cd.phenotype; 
	pNew->is_movable = true;
	pNew->is_out_of_domain = false;
	pNew->displacement.resize(3,0.0); // state? 
	
	pNew->assign_orientation();
	
	pNew->set_total_volume( pNew->phenotype.volume.total ); 
	
	return pNew; 
}

void Cell_Container_py::delete_cell( int index ){
    	
	PhysiCell::Cell* pDeleteMe = (*all_cells)[index]; 
	
	// release any attached cells (as of 1.7.2 release)
	pDeleteMe->remove_all_attached_cells(); 
	
	// released internalized substrates (as of 1.5.x releases)
	pDeleteMe->release_internalized_substrates(); 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	// shrink the vector
	(*all_cells).pop_back();	
	
	// deregister agent in from the agent container
	pDeleteMe->get_container()->remove_agent(pDeleteMe);
	// de-allocate (delete) the cell; 
	delete pDeleteMe; 


	return; 
    
}
void Cell_Container_py::delete_cell( PhysiCell::Cell* pDelete){
    
    delete_cell(pDelete->index);
    return; 
}

void Cell_Container_py::kill_cell(PhysiCell::Cell* pDelete){
    delete_cell(pDelete);
    return;
}

PhysiCell::Cell* Cell_Container_py::divide_cell(PhysiCell::Cell* pDivide){
    pDivide->remove_all_attached_cells(); 
	
	PhysiCell::Cell* child = create_cell();
	child->copy_data( pDivide );	
	child->copy_function_pointers(pDivide);
	child->parameters = pDivide->parameters;
	
	// evenly divide internalized substrates 
	// if these are not actively tracked, they are zero anyway 
	*pDivide->internalized_substrates *= 0.5; 
	*(child->internalized_substrates) = *pDivide->internalized_substrates ; 
	
	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );
	
	// randomly place the new agent close to me, accounting for orientation and 
	// polarity (if assigned)
		
	// May 30, 2020: 
	// Set cell_division_orientation = LegacyRandomOnUnitSphere to 
	// reproduce this code 
	/*
	double temp_angle = 6.28318530717959*UniformRandom();
	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
	
	double radius= phenotype.geometry.radius;
	std::vector<double> rand_vec (3, 0.0);
	
	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
	rand_vec[2]= cos( temp_phi );
	
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
	
	if( norm(rand_vec) < 1e-16 )
	{
		std::cout<<"************ERROR********************"<<std::endl;
	}
	normalize( &rand_vec ); 
	rand_vec *= radius; // multiply direction times the displacement 
	*/
	
	std::vector<double> rand_vec = PhysiCell::cell_division_orientation(); 
	rand_vec = rand_vec- pDivide->phenotype.geometry.polarity*(rand_vec[0]*pDivide->state.orientation[0]+ 
		rand_vec[1]*pDivide->state.orientation[1]+rand_vec[2]*pDivide->state.orientation[2])*pDivide->state.orientation;	
	rand_vec *= pDivide->phenotype.geometry.radius;

	child->assign_position(pDivide->position[0] + rand_vec[0],
						   pDivide->position[1] + rand_vec[1],
						   pDivide->position[2] + rand_vec[2]);
						 
	//change my position to keep the center of mass intact 
	// and then see if I need to update my voxel index
	static double negative_one_half = -0.5; 
	axpy( &pDivide->position, negative_one_half , rand_vec ); // position = position - 0.5*rand_vec; 

	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
	if( underlying_mesh.is_position_valid(pDivide->position[0], pDivide->position[1], pDivide->position[2]))
	{
		pDivide->is_out_of_domain = true;
		pDivide->is_active = false;
		pDivide->is_movable = false;
	}	
	 
	pDivide->update_voxel_in_container();
	pDivide->phenotype.volume.divide(); 
	child->phenotype.volume.divide();
	child->set_total_volume(child->phenotype.volume.total);
	pDivide->set_total_volume(pDivide->phenotype.volume.total);
	
	// child->set_phenotype( phenotype ); 
	child->phenotype = pDivide->phenotype; 
	
	return child;
}

void Cell_Container_py::save_all_cells_to_matlab( std::string filename ){}
void Cell_Container_py::delete_cell_original( int index ){
    //	std::cout << __FUNCTION__ << " " << (*all_cells)[index] 
//	<< " " << (*all_cells)[index]->type_name << std::endl; 
	
	// release any attached cells (as of 1.7.2 release)
	(*all_cells)[index]->remove_all_attached_cells(); 
	
	// released internalized substrates (as of 1.5.x releases)
	(*all_cells)[index]->release_internalized_substrates(); 
	
	// deregister agent in from the agent container
	(*all_cells)[index]->get_container()->remove_agent((*all_cells)[index]);
	// de-allocate (delete) the cell; 
	delete (*all_cells)[index]; 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	// shrink the vector
	(*all_cells).pop_back();	
	return; 
}




void save_PhysiCell_to_MultiCellDS_xml_pugi_py( std::string filename_base , 
												BioFVM_py::Microenvironment_py& M , 
												Cell_Container_py& CellCon, 
												double current_simulation_time)
{
    // start with a standard BioFVM save

    add_BioFVM_to_open_xml_pugi( BioFVM::biofvm_doc, filename_base , current_simulation_time , M ); 
    // now, add the PhysiCell data using the new container linked version

    add_PhysiCell_cells_to_open_xml_pugi_py( BioFVM::biofvm_doc , filename_base , M , CellCon); 
    // Lastly, save to the indicated filename 

    char filename[1024]; 
    sprintf( filename , "%s.xml" , filename_base.c_str() ); 
    BioFVM::biofvm_doc.save_file( filename );
	

    return; 
    
    
}


void add_PhysiCell_cells_to_open_xml_pugi_py( pugi::xml_document& xml_dom, std::string filename_base, Microenvironment& M, Cell_Container_py& CellCon  ){
	static double temp_zero = 0.0; 

	if( BioFVM::save_cell_data == false )
	{ return; }

	pugi::xml_node root = xml_dom.child("MultiCellDS") ; 
	pugi::xml_node node = root.child( "cellular_information" ); 
	root = node; 

	// Let's reduce memory allocations and sprintf calls. 
	// This reduces execution time by around 30%. (e.g., write time for 1,000,000 cells decreases from 
	// 45 to 30 seconds on an older machine. 
	static char* temp; 
	static bool initialized = false; 

	static char rate_chars [1024]; 
	static char volume_chars [1024]; 
	static char diffusion_chars [1024]; 
	if( !initialized )
	{ 
		temp = new char [1024]; 
		initialized = true; 
		
		sprintf( rate_chars, "1/%s" , M.time_units.c_str() ); 
		sprintf( volume_chars, "%s^3" , M.spatial_units.c_str() ); 
		sprintf( diffusion_chars , "%s^2/%s", M.spatial_units.c_str() , M.time_units.c_str() ); 
	}
	

	node = node.child( "cell_populations" ); 
	if( !node )
	{
		node = root.append_child( "cell_populations" ); 
	}
	root = node; // root = cell_populations 
	
	// if we are using the customized matlab data, do it here. 
	if( BioFVM::save_cells_as_custom_matlab == true || 1 == 1 )
	{
		node = node.child( "cell_population" ); 
		if( !node )
		{
			node = root.append_child( "cell_population" ); 
			pugi::xml_attribute attrib = node.append_attribute( "type" ); 
			attrib.set_value( "individual" ); 
		}
		
		if( !node.child( "custom" ) ) 
		{
			node.append_child( "custom" ); 
		}
		node = node.child( "custom" ); 
		
		// look for a node called simplified_data, with source = PhysiCell 
		
		pugi::xml_node node_temp = node.child( "simplified_data" ); 
		bool temp_search_done = false;
		while( !temp_search_done && node_temp )
		{
			if( node_temp )
			{
				pugi::xml_attribute attribute_temp = node_temp.attribute( "source" ); 
				if( attribute_temp )
				{
					if( strcmp( attribute_temp.value() , "PhysiCell" ) == 0 )
					{
						temp_search_done = true; 
					}
					else
					{
						node_temp = node_temp.next_sibling(); 
					}
				}
			}
			else
			{
				node_temp = (pugi::xml_node) NULL; 
			}
		}

		if( !node_temp )
		{
			node_temp = node.append_child( "simplified_data" ); 
			pugi::xml_attribute attrib = node_temp.append_attribute( "type" ); 
			attrib.set_value( "matlab" ) ; 
			
			attrib = node_temp.append_attribute( "source" ); 
			attrib.set_value("PhysiCell"); 

			int index = 0; 
			int size = 1; 
			
			pugi::xml_node node_temp1 = node_temp.append_child( "labels" ); 
			
			// ID,x,y,z,total volume
			node_temp1 = node_temp1.append_child( "label" ); 
			node_temp1.append_child( pugi::node_pcdata ).set_value( "ID" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 

			size = 3; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "position" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 

			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "total_volume" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
		
			// type, cycle model, current phase, elapsed time in phase, 
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "cell_type" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
		
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "cycle_model" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 					
			size = 1; 

			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "current_phase" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 			
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "elapsed_time_in_phase" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 			
			
			// nuclear volume, cytoplasmic volume, fluid fraction, calcified fraction, 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "nuclear_volume" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "cytoplasmic_volume" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "fluid_fraction" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 

			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "calcified_fraction" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			// orientation, polarity 			

			size = 3; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "orientation" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "polarity" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			// motility 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "migration_speed" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 

			size = 3; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "motility_vector" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "migration_bias" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 3; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "motility_bias_direction" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 
			
			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "persistence_time" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 			

			size = 1; 
			node_temp1 = node_temp1.append_child( "label" );
			node_temp1.append_child( pugi::node_pcdata ).set_value( "motility_reserved" ); 
			attrib = node_temp1.append_attribute( "index" ); 
			attrib.set_value( index ); 
			attrib = node_temp1.append_attribute( "size" ); 
			attrib.set_value( size ); 
			node_temp1 = node_temp1.parent(); 
			index += size; 			
			// custom variables 
std::cout<< "start";
			for( int i=0; i < (*CellCon.all_cells)[0]->custom_data.variables.size(); i++ )
			{
				size = 1; 
				char szTemp [1024]; 
std::cout<< " checkpoint 1";
				strcpy( szTemp, (*CellCon.all_cells)[0]->custom_data.variables[i].name.c_str() ); 
std::cout<< " checkpoint 2";
				node_temp1 = node_temp1.append_child( "label" );
				node_temp1.append_child( pugi::node_pcdata ).set_value( szTemp ); 
				attrib = node_temp1.append_attribute( "index" ); 
				attrib.set_value( index ); 
				attrib = node_temp1.append_attribute( "size" ); 
				attrib.set_value( size ); 
				node_temp1 = node_temp1.parent(); 
				index += size; 			
			}
			// custom vector variables 
			for( int i=0; i < (*CellCon.all_cells)[0]->custom_data.vector_variables.size(); i++ )
			{
				size = (*CellCon.all_cells)[0]->custom_data.vector_variables[i].value.size(); 
;				char szTemp [1024]; 
				strcpy( szTemp, (*CellCon.all_cells)[0]->custom_data.vector_variables[i].name.c_str() ); 
				node_temp1 = node_temp1.append_child( "label" );
				node_temp1.append_child( pugi::node_pcdata ).set_value( szTemp ); 
				attrib = node_temp1.append_attribute( "index" ); 
				attrib.set_value( index ); 
				attrib = node_temp1.append_attribute( "size" ); 
				attrib.set_value( size ); 
				node_temp1 = node_temp1.parent(); 
				index += size; 			
			}
std::cout<< "2";
			
		}
		std::cout<< "24";
		node = node_temp; 
		
		if( !node.child( "filename" ) )
		{
			node.append_child( "filename" ); 
		}
		node = node.child( "filename" ); 
		
		// next, filename 
		char filename [1024]; 
		sprintf( filename , "%s_cells_physicell.mat" , filename_base.c_str() ); 
		
		/* store filename without the relative pathing (if any) */ 
		char filename_without_pathing [1024];
		char* filename_start = strrchr( filename , '/' ); 
		if( filename_start == NULL )
		{ filename_start = filename; }
		else	
		{ filename_start++; } 
		strcpy( filename_without_pathing , filename_start );  
		
		if( !node.first_child() )
		{
			node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
		}
		else
		{
			node.first_child().set_value( filename_without_pathing ); // filename ); 
		}
		
		// next, create a matlab structure and save it!
		
		// order: ID,x,y,z,total volume, (same as BioFVM custom data, but instead of secretions ...)
		// type, cycle model, current phase, elapsed time in phase, 
		// nuclear volume, cytoplasmic volume, fluid fraction, calcified fraction, 
		// orientation, polarity 
		
		int number_of_data_entries = (*CellCon.all_cells).size(); 
		int size_of_each_datum = 1 + 3 + 1  // ID, x,y,z, total_volume 
			+1+1+1+1 // cycle information 
			+1+1+1+1 // volume information 
			+3+1 // orientation, polarity; 
			+1+3+1+3+1+1; // motility 
		// figure out size of custom data. for now, 
		// assume all the cells have teh same custom data as 
		// cell #0
		int custom_data_size = (*CellCon.all_cells)[0]->custom_data.variables.size();  
		for( int i=0; i < (*CellCon.all_cells)[0]->custom_data.vector_variables.size(); i++ )
		{
			custom_data_size += (*CellCon.all_cells)[0]->custom_data.vector_variables[i].value.size(); 
		}
		size_of_each_datum += custom_data_size; 
		

		FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "cells" );  
		if( fp == NULL )
		{ 
			std::cout << std::endl << "Error: Failed to open " << filename << " for MAT writing." << std::endl << std::endl; 
	
			std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
			<< "Check to make sure your save directory exists. " << std::endl << std::endl
			<< "I'm going to exit with a crash code of -1 now until " << std::endl 
			<< "you fix your directory. Sorry!" << std::endl << std::endl; 
			exit(-1); 
		} 
		PhysiCell::Cell* pCell; 
		
		// storing data as cols (each column is a cell)
		for( int i=0; i < number_of_data_entries ; i++ )
		{
			// ID, x,y,z, total_volume 
			double ID_temp = (double) (*CellCon.all_cells)[i]->ID;
			fwrite( (char*) &( ID_temp ) , sizeof(double) , 1 , fp ); 
			
			pCell = (*CellCon.all_cells)[i]; 

			fwrite( (char*) &( pCell->position[0] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->position[1] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->position[2] ) , sizeof(double) , 1 , fp ); 
			double dTemp = pCell->phenotype.volume.total; // get_total_volume();
			fwrite( (char*) &( dTemp ) , sizeof(double) , 1 , fp ); 
			
			// type, cycle model, current phase, elapsed time in phase, 
			dTemp = (double) pCell->type; 
			fwrite( (char*) &( dTemp ) , sizeof(double) , 1 , fp );  // cell type 
			
			dTemp = (double) pCell->phenotype.cycle.model().code; 
			fwrite( (char*) &( dTemp ) , sizeof(double) , 1 , fp ); // cycle model 
			
			dTemp = (double) pCell->phenotype.cycle.current_phase().code; 
			fwrite( (char*) &( dTemp ) , sizeof(double) , 1 , fp ); // current phase 
			
			// dTemp = pCell->phenotype.cycle.phases[pCell->phenotype.current_phase_index].elapsed_time; 
			fwrite( (char*) &( pCell->phenotype.cycle.data.elapsed_time_in_phase ) , sizeof(double) , 1 , fp ); // elapsed time in phase 
			
			// volume information
			// nuclear volume, cytoplasmic volume, fluid fraction, calcified fraction, 
			fwrite( (char*) &( pCell->phenotype.volume.nuclear ) , sizeof(double) , 1 , fp );  // nuclear volume 

			fwrite( (char*) &( pCell->phenotype.volume.cytoplasmic ) , sizeof(double) , 1 , fp );  // cytoplasmic volume 
			
			fwrite( (char*) &( pCell->phenotype.volume.fluid_fraction ) , sizeof(double) , 1 , fp );  // fluid fraction 

			fwrite( (char*) &( pCell->phenotype.volume.calcified_fraction ) , sizeof(double) , 1 , fp );  // calcified fraction 
			
			// orientation, polarity; 
			fwrite( (char*) &( pCell->state.orientation[0] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->state.orientation[1] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->state.orientation[2] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->phenotype.geometry.polarity ) , sizeof(double) , 1 , fp ); 
			

			
			// motility information 
			fwrite( (char*) &( pCell->phenotype.motility.migration_speed ) , sizeof(double) , 1 , fp ); // speed
			fwrite( (char*) &( pCell->phenotype.motility.motility_vector[0] ) , sizeof(double) , 1 , fp ); // velocity 
			fwrite( (char*) &( pCell->phenotype.motility.motility_vector[1] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->phenotype.motility.motility_vector[2] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->phenotype.motility.migration_bias ) , sizeof(double) , 1 , fp );  // bias (0 to 1)
			fwrite( (char*) &( pCell->phenotype.motility.migration_bias_direction[0] ) , sizeof(double) , 1 , fp ); // bias direction 
			fwrite( (char*) &( pCell->phenotype.motility.migration_bias_direction[1] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->phenotype.motility.migration_bias_direction[2] ) , sizeof(double) , 1 , fp ); 
			fwrite( (char*) &( pCell->phenotype.motility.persistence_time ) , sizeof(double) , 1 , fp ); // persistence 
			fwrite( (char*) &( temp_zero ) , sizeof(double) , 1 , fp ); // reserved for "time in this direction" 
			
			// custom variables 
			for( int j=0 ; j < pCell->custom_data.variables.size(); j++ )
			{
				fwrite( (char*) &( pCell->custom_data.variables[j].value ) , sizeof(double) , 1 , fp );  
			}
			
			// custom vector variables 
			for( int j=0 ; j < pCell->custom_data.vector_variables.size(); j++ )
			{
				for( int k=0; k < pCell->custom_data.vector_variables[j].value.size(); k++ )
				{
					fwrite( (char*) &( pCell->custom_data.vector_variables[j].value[k] ) , sizeof(double) , 1 , fp );  
				}
			}
			
		}

		fclose( fp ); 
		
		
		return; 
	}
	
	// if there's a list, clear it out 
	node = node.child( "cell_population" ); 
	if( node )
	{
		node = node.parent(); 
		node.remove_child( node.child( "cell_population" ) ); 
	}
	node = root.append_child( "cell_population" ); 
	pugi::xml_attribute attrib = node.append_attribute( "type" ); 
	attrib.set_value( "individual" ); 

	// now go through all cells 
	std::cout<< "3";
	root = node; 
	for( int i=0; i < (*CellCon.all_cells).size(); i++ )
	{
		node = node.append_child( "cell" ); 
		attrib = node.append_attribute( "ID" ); 
		attrib.set_value(  (*CellCon.all_cells)[i]->ID ); 
		
		node = node.append_child( "phenotype_dataset" ); 
		node = node.append_child( "phenotype" ); // add a type? 
		
		// add all the transport information 
		node = node.append_child( "transport_processes" ); 
		
		// add variables and their source/sink/saturation values (per-cell basis)
		for( int j=0; j < M.number_of_densities() ; j++ ) 
		{
			node = node.append_child( "variable" ); 
			attrib = node.append_attribute( "name" ); 
			attrib.set_value( M.density_names[j].c_str() ); 
			// ChEBI would go here later 
			attrib = node.append_attribute( "ID" );
			attrib.set_value( j ); 
			
			node = node.append_child( "export_rate" ); 
			attrib = node.append_attribute( "units" ); 
			attrib.set_value( rate_chars ); 
			// sprintf( temp , "%f" , all_basic_agents[i]->get_total_volume() * (*all_basic_agents[i]->secretion_rates)[j] ); 
			sprintf( temp , "%f" , (*CellCon.all_cells)[i]->phenotype.volume.total * (*(*CellCon.all_cells)[i]->secretion_rates)[j] ); 
			node.append_child( pugi::node_pcdata ).set_value( temp ); 
			node = node.parent( ); 
			
			node = node.append_child( "import_rate" ); 
			attrib = node.append_attribute( "units" ); 
			attrib.set_value( rate_chars ); 
			sprintf( temp,  "%f" , (*CellCon.all_cells)[i]->phenotype.volume.total * (*(*CellCon.all_cells)[i]->uptake_rates)[j] ); 
			node.append_child( pugi::node_pcdata ).set_value( temp ); 
			node = node.parent(); 
			
			node = node.append_child( "saturation_density" ); 
			attrib = node.append_attribute( "units" ); 
			attrib.set_value( M.density_units[j].c_str() ); 
			sprintf( temp, "%f" , (*(*CellCon.all_cells)[i]->saturation_densities)[j] ); 
			node.append_child( pugi::node_pcdata ).set_value( temp ); 
			node = node.parent(); 
			
			node = node.parent(); // back up to transport processes 
			std::cout<< "5";
		}
		
		node = node.parent(); // back up to phenotype 

		// add size information 
		node = node.append_child( "geometrical_properties" ); 
		node = node.append_child("volumes");
		node = node.append_child("total_volume"); 
		attrib = node.append_attribute("units"); 
		attrib.set_value( volume_chars ); 
		sprintf( temp,  "%f" , (*CellCon.all_cells)[i]->phenotype.volume.total ); 
		node.append_child( pugi::node_pcdata ).set_value( temp ); 		
		node = node.parent(); 
		node = node.parent(); 
		
		node = node.parent(); // back up to geometrical_properties 
		
		node = node.parent(); // back up to phenotype 
		
		node = node.parent(); // back up to phenotype_dataset 
		
		// add position information 
		node = node.append_child( "state"); 
		node = node.append_child( "position" ); 
		attrib = node.append_attribute( "units" ); 
		attrib.set_value( M.spatial_units.c_str() ); 
		
		// vector3_to_list( all_basic_agents[i]->position , temp , ' ');
		sprintf( temp , "%.7e %.7e %.7e" , (*CellCon.all_cells)[i]->position[0], (*CellCon.all_cells)[i]->position[1], (*CellCon.all_cells)[i]->position[2] ); 
		node.append_child( pugi::node_pcdata ).set_value( temp ); 		
		
		node = root;
		std::cout<< "6";
	}
	
	return; 
}
  
    
    
}//end namespace PhysiCellCore_py
