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
    
    #setting the up the microenvironment
    mcrenv = setup_microenvironment(mcrenv,settings)
    voxelsize = 30
    
    #setting up the cell container
    cellcontainer = core.Cell_Container(microenvironment = mcrenv, voxel_size = voxelsize)
    current_time = 0
    max_time = 5.0
    
    diffusion_dt = 1.0
    
    mcrenv.simulate_diffusion_decay(diffusion_dt)
    #Update all the cells
    cellcontainer.update_all_cells(current_time, diffusion_dt, diffusion_dt, diffusion_dt)
    
    
    while(current_time < (max_time + diffusion_dt)):
         #Simulate Diffusion in the microenvironment
         mcrenv.simulate_diffusion_decay(diffusion_dt)
         #Update all the cells
         cellcontainer.update_all_cells(current_time, diffusion_dt, diffusion_dt, diffusion_dt)
         #Update the current time
         core.save_MultiCellDS(filename = "test"+str(current_time), microenvironment = mcrenv, cell_container = cellcontainer, current_simulation_time = current_time)
         current_time += diffusion_dt