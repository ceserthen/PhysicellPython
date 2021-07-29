# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:44:49 2021

@author: belch
"""
# Start of test module for heterogeneity
import os
os.environ["OMP_NUM_THREADS"] = "4"
import sys
sys.path.append('../../lib')
from physicell import BioFVM
from physicell import core

#Initialize the microenvironment with parameters from the xml file
mcrenv = BioFVM.Microenvironment(xml_config_file ="PhysiCell_settings.xml")

mcrenv.initialize_microenvironment()
mcrenv.name = "synthetic tissue"
mcrenv.find_density_index("cargo signal")