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
from math import fabs
from modules import multicellds
from xml_parse import PhysicellSettings
from visualization import visualize_microenvironment
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
#def visualize_microenvironment():
    
#biorobots main
if __name__ == "__main__":
    xml_config_file ="PhysiCell_settings.xml"
    mcrenv = BioFVM.Microenvironment(xml_config_file)
    settings = PhysicellSettings(xml_config_file)
    mcrenv = setup_microenvironment(mcrenv,settings)
    voxelsize = 30.0
    cellcontainer = core.Cell_Container(microenvironment = mcrenv, voxel_size = voxelsize)
    #Following line crashes
    #    cellcontainer.create_cell()
    current_time = 0.0
#    max_time = 30000.0
#    diffusion_dt = 1.0
    #get max_time diffusion_dt from file
    max_time = settings.save_settings.max_time
    diffusion_dt= settings.save_settings.dt_diffusion
    mcrenv.simulate_diffusion_decay(diffusion_dt)

    cellcontainer.update_all_cells(current_time, diffusion_dt, diffusion_dt, diffusion_dt)
    idx=0
    next_full_save_time = 0.0
    while(current_time < (max_time + diffusion_dt)):
        
         #Simulate Diffusion in the microenvironment
        #update microenvironment
        mcrenv.simulate_diffusion_decay(diffusion_dt)
        #update cells
        cellcontainer.update_all_cells(current_time, diffusion_dt, diffusion_dt, diffusion_dt)
        #Update the current time
        current_time += diffusion_dt
        idx+=1

        #save if its time        
        if((current_time - next_full_save_time) < 0.01*diffusion_dt ):
            next_full_save_time += float(settings.save_settings.interval)
            print(settings.save_settings.interval)
#            core.save_MultiCellDS(filename = "output0000000"+str(idx), microenvironment = mcrenv, cell_container = cellcontainer, current_simulation_time = current_time)
        
        mcrenv.display_info()
    visualize_microenvironment(4)
    visualize_microenvironment(5)