#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 00:42:14 2021

@author: bscuser
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tdiniaco
"""

import xml.etree.ElementTree as ET
class SaveSettings(object):
    def __init__(self,tree):
#        tree = ET.parse(settings_xml_fname)
        root = tree.getroot()
        
        #overall node
        overall_node = root.findall("overall")[0]
        node = overall_node.findall("max_time")[0]
        self._max_time = int(float(node.text))
        self._max_time_units = node.attrib['units']
        #save node
        save_node = root.findall("save")[0]
        node = save_node.findall("full_data")[0]
        node = node.findall("interval")[0]
        self._interval = int(float(node.text))
        self._interval_units = node.attrib['units']
        node = overall_node.findall("dt_diffusion")[0]
        self.dt_diffusion = float(node.text)
        node = overall_node.findall("dt_mechanics")[0]
        self.dt_mechanics = float(node.text)
        node = overall_node.findall("dt_phenotype")[0]
        self.dt_phenotype = float(node.text)
    @property
    def max_time(self):
        return self._max_time
    @property
    def max_time_units(self):
        return self._max_time_units
    @property
    def interval(self):
        return self._interval
    @property
    def interval_units(self):
        return self._interval_units

class UserParameterSettings(object):      
    def __init__(self,tree):
        root = tree.getroot()
        parameter_node = root.findall("user_parameters")[0]
#        print(root.tag)
#        print(parameter_node.tag)
        self._parameter_dict=self.set_parameters(parameter_node)
        #overall node
    def set_parameters(self,root):
        params = {}
        for node in root:
            param= Parameters(node.attrib['type'],node.attrib['units'],node.text)
            params[node.tag]= param
        return params
    @property
    def parameter_dict(self):
        return self._parameter_dict
    def parameter_dict_value(self,name):
        return self.parameter_dict[name].parameter_value      
class Parameters(object):
    def __init__(self,ptype,punits,pvalue):
        self._parameter_type=ptype
        self._parameter_units=punits
        self._parameter_value=pvalue
    @property
    def parameter_type(self):
        return self._parameter_type   
    @property
    def parameter_units(self):
        return self._parameter_units
    @property
    def parameter_value(self):
        if self._parameter_type=="int":
            self._parameter_value = int(self._parameter_value)
        if self._parameter_type=="double":
            self._parameter_value = float(self._parameter_value)
        if self._parameter_type=="string":
            self._parameter_valueint = str(self._parameter_value)
        return self._parameter_value           
        #
class PhysicellSettings(object):
    def __init__(self,settings_xml_fname='PhysiCell_settings.xml'):
        self._tree = ET.parse(settings_xml_fname)
        root = self._tree.getroot()
        self.user_parameters=UserParameterSettings(self._tree)
        self.save_settings = SaveSettings(self._tree)
#    def get_param_value(self):
        
        