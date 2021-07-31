# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 23:56:31 2021

Legacy variant of pyPhysicell originally by @rheiland
https://github.com/rheiland/PyPhysiCel
@author: Donald Belcher
"""

import os
os.environ["OMP_NUM_THREADS"] = "4"
import sys
sys.path.append('../../lib')
from physicell import BioFVM
from physicell import core
from physicell import biorobots
from modules import multicellds
from xml_parse import PhysicellSettings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection
#class biorobots:

time_delay = 0.1

#-----------------------------------------------------
def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83 
    Make a scatter plot of circles. 
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_)
               for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    return collection


fig = plt.figure(figsize=(7,7))
ax = fig.gca()

# --- TODO! Delete model, e.g. all_cells->size()  --> 0

# retval = load_PhysiCell_config_file("config/cancerbots.xml")
retval = core.load_PhysiCell_config_file("biorobots.xml")

BioFVM.initialize_microenvironment_leg()   # in setup_microenvironment()
# double mechanics_voxel_size = 30; 
# this creates  std::vector<Cell*> *all_cells;
mechanics_voxel_size = 30
# --- original method took pointer to menv!
# Cell_Container* create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size )
# --- 
# Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
#create_cell_container_for_microenvironment( menv, mechanics_voxel_size )
core.create_cell_container2( mechanics_voxel_size )

biorobots.create_cell_types_biorobots()
# setup_tissue1()
biorobots.setup_tissue_biorobots()

settings = pcore.get_physicell_settings()
print('settings.omp_num_threads=',settings.omp_num_threads)
print('settings.max_time=',settings.max_time)

pc_globals = pcore.get_physicell_globals()

# rf. PhysiCell_constants.h
diffusion_dt = 0.01
mechanics_dt = 0.1
phenotype_dt = 6.0

#tdel = settings.max_time/4.0

# display_citations(); 
## set the performance timers 
# BioFVM::RUNTIME_TIC();
# BioFVM::TIC();

# for mytime in np.arange(0, settings.max_time, tdel):
max_steps = settings.max_time / diffusion_dt
for istep in range(0, 4):
    # print('time (for loop)=',mytime)
    print('current_time=',pc_globals.current_time)

    # update the microenvironment
	# microenvironment.simulate_diffusion_decay( diffusion_dt );
    menv.simulate_diffusion_decay( diffusion_dt )

	# run PhysiCell 
	# ((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
#void Cell_Container::update_all_cells(double t)
#{
#	update_all_cells(t, phenotype_dt, mechanics_dt , diffusion_dt );
    # update_all_cells(mytime)
    update_all_cells(pc_globals.current_time)
    print('# cells= ',get_num_cells())

    # sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
    # SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
			
	# PhysiCell_globals.current_time += diffusion_dt
    pc_globals.current_time += diffusion_dt

    pos = get_cells_pos2D()
    plt.cla()
    # title_str += " (" + str(num_cells) + " agents)"
    title_str = " (" + str(len(pos)) + " agents)"
    # plt.title(title_str)
    # plt.xlim(axes_min,axes_max)
    # plt.ylim(axes_min,axes_max)
    #  plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs)
    #  plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs, alpha=0.5, edgecolor='black')
    #  plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs, alpha=1.0, edgecolor='black')
    #  circles(xvals,yvals, s=rvals, c=rgbs, alpha=1.0, edgecolor='black')
    #  circles(xvals,yvals, s=rvals)
    #  circles(xvals,yvals, s=rvals, c=rgbs)

    # circles(xvals,yvals, s=rvals, color=rgbs)
    pos = get_cells_x()
    xvals = np.asarray(pos)
    pos = get_cells_y()
    yvals = np.asarray(pos)
#    circles(xvals,yvals, s=2.5)

    #plt.xlim(0,2000)  # TODO - get these values from width,height in .svg at top
    #plt.ylim(0,2000)

    plt.pause(time_delay)

pos = get_cells_x()
xvals = np.asarray(pos)
pos = get_cells_y()
yvals = np.asarray(pos)

pos = get_cells_types()
ctype = np.asarray(pos)

#rgb_tuple = mplc.to_rgb(mplc.cnames[s])
red_tuple = mplc.to_rgb(mplc.cnames['red'])   # (1.0, 0.0, 0.0)
blue_tuple = mplc.to_rgb(mplc.cnames['blue'])

rgb_list = []
#rgb_list.append(rgb)
rgb_list.append(red_tuple)
cvals = [blue_tuple for x in pos]
for idx in range(len(cvals)):
    if pos[idx] == 0:
        cvals[idx] = red_tuple

rgbs = np.asarray(cvals)

circles(xvals,yvals, s=7.0, color=rgbs)
title_str = " biorobots with pyphysicell: " + str(len(pos)) + " agents"
plt.title(title_str)
#pos = get_cells_pos2D()
#a = np.asarray(pos)  # len(a)