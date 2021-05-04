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
#ifndef __BioFVM_py_h__
#define __BioFVM_py_h__


#include "../../../src/BioFVM/BioFVM.h"
#include "../../../src/BioFVM/pugixml.hpp"
#include "../../../src/modules/PhysiCell_settings.h"
#include "../../../src/modules/PhysiCell_pugixml.h"

namespace BioFVM_py{
    
  class Microenvironment_py : public BioFVM::Microenvironment {
   
public:
    BioFVM::Microenvironment_Options microenvironment_options;
    
    Microenvironment_py();
    Microenvironment_py( BioFVM::Microenvironment_Options const microenvironment_options_in );
    Microenvironment_py(std::string filename);
    Microenvironment_py( pugi::xml_node root_node );
    
    bool setup_microenvironment_from_XML_node( pugi::xml_node root_node ); 
    void setup_microenvironment_from_file( std::string filename );
    
    //Microenvironment_initializer
    void initialize_microenvironment( void );
    void display_info (void);
    
    //overwriting for internal handling of microenvironemtn options
    void display_information( std::ostream& os );
    void resize_densities( int new_size );
    void add_density( void );
    void add_density( std::string name , std::string units );
    void add_density( std::string name , std::string units, double diffusion_constant, double decay_rate );
    void add_density_helper(void);
    void set_density( int index , std::string name , std::string units );
    void set_density( int index , std::string name , std::string units , double diffusion_constant , double decay_rate );
    
};  

    
}

#endif
