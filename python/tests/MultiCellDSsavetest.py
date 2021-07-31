# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:44:49 2021

@author: belch
"""
# Start of test module for heterogeneity
import os
os.environ["OMP_NUM_THREADS"] = "4"
import sys
import scipy.io
sys.path.append('../../lib')
from physicell import BioFVM
from physicell import core

#Initialize the microenvironment with parameters from the xml file
mcrenv = BioFVM.Microenvironment(xml_config_file ="PhysiCell_settings.xml")

mcrenv.initialize_microenvironment()
mcrenv.name = "synthetic tissue"
mcrenv.find_density_index("cargo signal")

voxelsize = 30.0
cellcontainer = core.Cell_Container(microenvironment = mcrenv, voxel_size = voxelsize)

for i in [1,2,3]:
    mcrenv.simulate_diffusion_decay(1)
    core.save_MultiCellDS(filename = "test"+str(i), microenvironment = mcrenv, cell_container = cellcontainer, current_simulation_time = i)
