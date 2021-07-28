# Start of test module for heterogeneity
import os
os.environ["OMP_NUM_THREADS"] = "4"
import sys
sys.path.append('../../lib')
from physicell import BioFVM
from physicell import core

#Initialize the microenvironment with parameters from the xml file
mcrenv = BioFVM.Microenvironment(xml_config_file ="defaultHeterogeneity_settings.xml")

#Setting the voxel size
voxelsize = 30.0

#Initializing the Cell Container using the Microenvironment and Voxel Size
cellcontainer = core.Cell_Container(microenvironment = mcrenv, voxel_size = voxelsize)
current_time = 0
max_time = 30000.0

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
     current_time += diffusion_dt