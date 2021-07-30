#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 21:00:45 2021

@author: bscuser
"""

# Start of test module for biorobots
import os
os.environ["OMP_NUM_THREADS"] = "4"
import sys
sys.path.append('../../lib')
from physicell import BioFVM
from physicell import core
from modules import multicellds
from xml_parse import PhysicellSettings
#class biorobots:
def setup_microenvironment(microenv,settings):
    
    microenv.initialize_microenvironment()
    microenv.name = "synthetic tissue"
    print(mcrenv.name)
    cargo_index = microenv.find_density_index("cargo signal")
    director_index = microenv.find_density_index("director signal")
    microenv.diffusion_coefficients[cargo_index] = settings.user_parameters.parameter_dict_value("cargo_signal_D")#		parameters.doubles("director_signal_D");  
    microenv.decay_rates[cargo_index] = settings.user_parameters.parameter_dict_value("cargo_signal_decay")

	microenv.diffusion_coefficients[director_index] = settings.user_parameters.parameter_dict_value("director_signal_D")
    microenv.decay_rates[director_index] = settings.user_parameters.parameter_dict_value("director_signal_decay")#		parameters.doubles("director_signal_D");  

    return microenv
if __name__ == "__main__":
    xml_config_file ="PhysiCell_settings.xml"
    mcrenv = BioFVM.Microenvironment(xml_config_file)
    settings = PhysicellSettings(xml_config_file)
    
    mcrenv = setup_microenvironment(mcrenv,settings)