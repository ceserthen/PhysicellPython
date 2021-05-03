
#include "../../../pybind11/include/pybind11/pybind11.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
 #include <pybind11/stl.h>
 #include <pybind11/stl_bind.h>


#include <iostream>
#include <vector> 


#include "../../../src/BioFVM/BioFVM.h"
#include "../../../src/modules/PhysiCell_settings.h"

namespace py = pybind11;


PYBIND11_MAKE_OPAQUE(std::vector<std::string, std::allocator<std::string>>);
using StringList = std::vector<std::string, std::allocator<std::string>>;

PYBIND11_MAKE_OPAQUE(std::vector<double, std::allocator<double>>);
using DoubleList = std::vector<double, std::allocator<double>>;
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);

namespace BioFVM{
class Microenvironment_py : public Microenvironment {
private:
    
public:
    Microenvironment_Options microenvironment_options;
    
    Microenvironment_py();
    Microenvironment_py( Microenvironment_Options const &microenvironment_options_in );
    Microenvironment_py( std::string name, std::string filename);
    Microenvironment_py( pugi::xml_node root_node );
    
    bool setup_microenvironment_from_XML_node( pugi::xml_node root_node );
    bool setup_microenvironment_from_XML( std::string filename );
    
    //overwriting for internal handling of microenvironemtn options
    //TODO: Find these dependecies to see if I need to overwrite any other internal functions
    void display_information( std::ostream& os ); 
    void resize_densities( int new_size );
    void add_density( void );
    void add_density( std::string name , std::string units );
    void add_density( std::string name , std::string units, double diffusion_constant, double decay_rate );
    void set_density( int index , std::string name , std::string units );
    void set_density( int index , std::string name , std::string units , double diffusion_constant , double decay_rate );
    
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
    diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 

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

    default_microenvironment_options.Dirichlet_all.assign( 1 , true ); 
    default_microenvironment_options.Dirichlet_xmin.assign( 1 , false ); 
    default_microenvironment_options.Dirichlet_xmax.assign( 1 , false ); 
    default_microenvironment_options.Dirichlet_ymin.assign( 1 , false ); 
    default_microenvironment_options.Dirichlet_ymax.assign( 1 , false ); 
    default_microenvironment_options.Dirichlet_zmin.assign( 1 , false ); 
    default_microenvironment_options.Dirichlet_zmax.assign( 1 , false ); 

    default_microenvironment_options.Dirichlet_xmin_values.assign( 1 , 1.0 ); 
    default_microenvironment_options.Dirichlet_xmax_values.assign( 1 , 1.0 ); 
    default_microenvironment_options.Dirichlet_ymin_values.assign( 1 , 1.0 ); 
    default_microenvironment_options.Dirichlet_ymax_values.assign( 1 , 1.0 ); 
    default_microenvironment_options.Dirichlet_zmin_values.assign( 1 , 1.0 ); 
    default_microenvironment_options.Dirichlet_zmax_values.assign( 1 , 1.0 ); 

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
    this->name=name;
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
    
    //simulation methods
    .def("simulate_diffusion_decay", static_cast<void (BioFVM::Microenvironment::*)(double)>(&BioFVM::Microenvironment::simulate_diffusion_decay), "advance the diffusion-decay solver by dt time");
}


