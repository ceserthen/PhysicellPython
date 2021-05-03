
#include "../../../pybind11/include/pybind11/pybind11.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#include <iostream>
#include <vector> 

#include "../../../src/BioFVM/BioFVM.h"
#include "../../../src/BioFVM/pugixml.hpp"
#include "../../../src/modules/PhysiCell_settings.h"
#include "../../../src/modules/PhysiCell_pugixml.h"


namespace py = pybind11;



namespace BioFVM_py
{

class Microenvironment_py : public BioFVM::Microenvironment {
   
public:
    BioFVM::Microenvironment_Options microenvironment_options;
    
    Microenvironment_py();
    Microenvironment_py( BioFVM::Microenvironment_Options const microenvironment_options_in );
    Microenvironment_py( std::string name, std::string filename);
    Microenvironment_py( pugi::xml_node root_node );
    
    bool setup_microenvironment_from_XML_node( pugi::xml_node root_node ); //TODO
    bool setup_microenvironment_from_XML( std::string filename );
    
    //overwriting for internal handling of microenvironemtn options
    //TODO: Find these dependecies to see if I need to overwrite any other internal functions
    void display_information( std::ostream& os ); //TODO
    void resize_densities( int new_size );//TODO
    void add_density( void );//TODO
    void add_density( std::string name , std::string units );//TODO
    void add_density( std::string name , std::string units, double diffusion_constant, double decay_rate );//TODO
    void set_density( int index , std::string name , std::string units );//TODO
    void set_density( int index , std::string name , std::string units , double diffusion_constant , double decay_rate );//TODO
    
};



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

    BioFVM::default_microenvironment_options.Dirichlet_all.assign( 1 , true ); 
    BioFVM::default_microenvironment_options.Dirichlet_xmin.assign( 1 , false ); 
    BioFVM::default_microenvironment_options.Dirichlet_xmax.assign( 1 , false ); 
    BioFVM::default_microenvironment_options.Dirichlet_ymin.assign( 1 , false ); 
    BioFVM::default_microenvironment_options.Dirichlet_ymax.assign( 1 , false ); 
    BioFVM::default_microenvironment_options.Dirichlet_zmin.assign( 1 , false ); 
    BioFVM::default_microenvironment_options.Dirichlet_zmax.assign( 1 , false ); 

    BioFVM::default_microenvironment_options.Dirichlet_xmin_values.assign( 1 , 1.0 ); 
    BioFVM::default_microenvironment_options.Dirichlet_xmax_values.assign( 1 , 1.0 ); 
    BioFVM::default_microenvironment_options.Dirichlet_ymin_values.assign( 1 , 1.0 ); 
    BioFVM::default_microenvironment_options.Dirichlet_ymax_values.assign( 1 , 1.0 ); 
    BioFVM::default_microenvironment_options.Dirichlet_zmin_values.assign( 1 , 1.0 ); 
    BioFVM::default_microenvironment_options.Dirichlet_zmax_values.assign( 1 , 1.0 ); 

    microenvironment_options = default_microenvironment_options;
    return;
};

Microenvironment_py::Microenvironment_py(BioFVM::Microenvironment_Options microenvironment_options_in)
{
    Microenvironment_py();
    microenvironment_options = microenvironment_options_in;
    return;
};
 

Microenvironment_py::Microenvironment_py( std::string name, std::string filename)
{

    Microenvironment_py(root_node);
    setup_microenvironment_from_XML(filename);
    this -> name = name;
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
		
		// now, figure out if individual boundaries are set 
/*		
		if( node1.attribute("boundaries") )
		{
			std::string option_string = node1.attribute("boundaries").value(); 
			Dirichlet_all.push_back(false); 

			if( strstr( option_string.c_str() , "xmin" ) )
			{ Dirichlet_xmin.push_back( true ); }
			else
			{ Dirichlet_xmin.push_back( false ); }
		
			if( strstr( option_string.c_str() , "xmax" ) )
			{ Dirichlet_xmax.push_back( true ); }
			else
			{ Dirichlet_xmax.push_back( false ); }
		
			if( strstr( option_string.c_str() , "ymin" ) )
			{ Dirichlet_ymin.push_back( true ); }
			else
			{ Dirichlet_ymin.push_back( false ); }
		
			if( strstr( option_string.c_str() , "ymax" ) )
			{ Dirichlet_ymax.push_back( true ); }
			else
			{ Dirichlet_ymax.push_back( false ); }
		
			if( strstr( option_string.c_str() , "zmin" ) )
			{ Dirichlet_zmin.push_back( true ); }
			else
			{ Dirichlet_zmin.push_back( false ); }

			if( strstr( option_string.c_str() , "zmax" ) )
			{ Dirichlet_zmax.push_back( true ); }
			else
			{ Dirichlet_zmax.push_back( false ); }
		}
		else
		{	
			Dirichlet_all.push_back(true); 
		}
*/		
		
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
	
	return true;  
};

bool Microenvironment_py::setup_microenvironment_from_XML( std::string filename )
{
    pugi::xml_document physicell_config_doc;
    pugi::xml_node physicell_config_root; 
    
    std::cout << "Using config file for microenvironment:" << filename << " ... " << std::endl ; 
    pugi::xml_parse_result result = physicell_config_doc.load_file( filename.c_str()  );

    if( result.status != pugi::xml_parse_status::status_ok )
    {
        std::cout << "Error loading " << filename << "!" << std::endl; 
        return false;
    }

    physicell_config_root = physicell_config_doc.child("PhysiCell_settings");
    setup_microenvironment_from_XML_node(physicell_config_root);
    return true;
};

}; //end namespace BioFVM

PYBIND11_MAKE_OPAQUE(std::vector<std::string, std::allocator<std::string>>);
using StringList = std::vector<std::string, std::allocator<std::string>>;

PYBIND11_MAKE_OPAQUE(std::vector<double, std::allocator<double>>);
using DoubleList = std::vector<double, std::allocator<double>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);


PYBIND11_MODULE(biofvmpy, m) 
{
    
    py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::bind_vector<std::vector<std::string>>(m, "VectorString");
    
    //TODO: Not currently working DAB 05-03-21
    //py::bind_vector<std::vector<std::vector<double>>>(m, "VectorVectorDouble");
    
    
    
    py::class_<BioFVM::Microenvironment>(m, "Microenvironment")
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
    
    .def("resize_voxels", static_cast<void (BioFVM::Microenvironment::*)(int)>(&BioFVM::Microenvironment::resize_voxels), "resize_voxels: enter the new number of voxels")
    
    .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(int, int, int)>(&BioFVM::Microenvironment::resize_space), "resize_space")
    .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(double, double, double, double, double, double, int, int, int)>(&BioFVM::Microenvironment::resize_space), "resize_space")
    .def("resize_space", static_cast<void (BioFVM::Microenvironment::*)(double,double,double,double,double,double,double,double,double)> (&BioFVM::Microenvironment::resize_space), "resize_space")
    
    .def("resize_space_uniform", static_cast<void (BioFVM::Microenvironment::*)(double,double,double,double,double,double,double)> (&BioFVM::Microenvironment::resize_space_uniform), "resize_space_uniform")
    
    //simulation methods
    .def("simulate_diffusion_decay", static_cast<void (BioFVM::Microenvironment::*)(double)>(&BioFVM::Microenvironment::simulate_diffusion_decay), "advance the diffusion-decay solver by dt time");
    
    
    py::class_<BioFVM_py::Microenvironment_py, BioFVM::Microenvironment>(m, "Microenvironmentpy")
    //constructors    
        .def(py::init<>())
        .def(py::init<BioFVM::Microenvironment_Options>())
        .def(py::init<std::string, std::string>())
        .def(py::init<pugi::xml_node>())
    
    //XML readers
        .def("setup_microenvironment_from_XML_node", static_cast<bool (BioFVM_py::Microenvironment_py::*)(pugi::xml_node)>(&BioFVM_py::Microenvironment_py::setup_microenvironment_from_XML_node), "Load microenvironment settings from an XML node ");
        
        //.def("setup_microenvironment_from_XML", static_cast<bool (BioFVM_py::Microenvironment_py::*)(std::string)>(&BioFVM_py::Microenvironment_py::setup_microenvironment_from_XML, "Load microenvironment settings from an XML File ");
}

// public:
//     BioFVM::Microenvironment_Options microenvironment_options;
//     
//     Microenvironment_py();
//     Microenvironment_py( BioFVM::Microenvironment_Options const microenvironment_options_in );
//     Microenvironment_py( std::string name, std::string filename);
//     Microenvironment_py( pugi::xml_node root_node );
//     
//     bool setup_microenvironment_from_XML_node( pugi::xml_node root_node ); //TODO
//     bool setup_microenvironment_from_XML( std::string filename );
//     
//     //overwriting for internal handling of microenvironemtn options
//     //TODO: Find these dependecies to see if I need to overwrite any other internal functions
//     void display_information( std::ostream& os ); //TODO
//     void resize_densities( int new_size );//TODO
//     void add_density( void );//TODO
//     void add_density( std::string name , std::string units );//TODO
//     void add_density( std::string name , std::string units, double diffusion_constant, double decay_rate );//TODO
//     void set_density( int index , std::string name , std::string units );//TODO
//     void set_density( int index , std::string name , std::string units , double diffusion_constant , double decay_rate );//TODO
//     
// };
