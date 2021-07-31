#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 12:09:30 2021

@author: bscuser
"""

import sys,pathlib
import xml.etree.ElementTree as ET
import math
import scipy.io
import matplotlib
import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt  # NB! do this AFTER the TkAgg line below!
#import matplotlib.colors as mplc
import numpy as np
from numpy.random import randn
def plot_substrate():
    global params
#    global current_idx, axes_max, cbar
    # select whichever substrate index you want, e.g., for one model:
    # 4=tumor cells field, 5=blood vessel density, 6=growth substrate
#    field_index = 4
#    field_index = 5
#    xml_file = "output%08d.xml" % current_idx
#    current_idx=0
    xml_file = "output00000001.xml"
    tree = ET.parse(xml_file)
    root = tree.getroot()
    print('time=' + root.find(".//current_time").text)
#     print('debug> plot_substrate: idx=',FileId)
#    fname = "output%08d_microenvironment0.mat" % FileId
    fname = "output%08d_microenvironment0.mat" % params["current_idx"]
    output_dir_str = '.'
    fullname = output_dir_str + "/" + fname
    if not pathlib.Path(fullname).is_file():
        print("file not found",fullname)
        return
    info_dict = {}
    scipy.io.loadmat(fullname, info_dict)
    M = info_dict['multiscale_microenvironment']
#     global_field_index = int(mcds_field.value)
    print('plot_substrate: field_index=',params['field_index'])
    f = M[params['field_index'],:]   #
    #plt.clf()
    #my_plot = plt.imshow(f.reshape(400,400), cmap='jet', extent=[0,20, 0,20])

#    fig = plt.figure(figsize=(7,7))
#    fig = plt.figure(figsize=(7,5.8))
    N = int(math.sqrt(len(M[0,:])))
    grid2D = M[0,:].reshape(N,N)
    xvec = grid2D[0,:]
    #xvec.size
    #xvec.shape
    num_contours = 30
    num_contours = 10
#    my_plot = plt.contourf(xvec,xvec,M[field_index,:].reshape(100,100), num_contours, cmap='viridis') #'viridis'
    my_plot = plt.contourf(xvec,xvec,M[params['field_index'],:].reshape(N,N), num_contours, cmap='viridis') #'viridis'
#    cbar.remove()
    if params["cbar"] == None:
      params["cbar"] = plt.colorbar(my_plot)
    else:
      params["cbar"] = plt.colorbar(my_plot, cax=params["cbar"].ax)
    axes_min = 0
    axes_min = -2000
    axes_max = 2000
#    plt.xlim(axes_min,axes_max)
#    plt.ylim(axes_min,axes_max)
#     plt.title(fname)
#     ax.set_title(fname)
#    plt.axis('equal')
#
# #    plt.show()
    plt.pause(params['time_delay'])
def press(event):
#  global current_idx, step_value
    global params
#    print('press', event.key)
    sys.stdout.flush()
    if event.key == 'escape':
        sys.exit(1)
    elif event.key == 'h':  # help
        print('esc: quit')
        print('right arrow: increment by step_value')
        print('left arrow:  decrement by step_value')
        print('up arrow:   increment step_value by 1')
        print('down arrow: decrement step_value by 1')
        print('0: reset to 0th frame')
        print('h: help')
    elif event.key == 'left':  # left arrow key
#    print('go backwards')
#    fig.canvas.draw()
        params["current_idx"] -= params["step_value"]
        if (params["current_idx"] < 0):
            params["current_idx"] = 0
        plot_substrate()
    elif event.key == 'right':  # right arrow key
#        print('go forwards')
#        fig.canvas.draw()
        params["current_idx"] += params["step_value"]
        plot_substrate()
    elif event.key == 'up':  # up arrow key
        params["step_value"] += 1
        print('step_value=',params["step_value"])
    elif event.key == 'down':  # down arrow key
        params["step_value"] -= 1
        if (params["step_value"] <= 0):
            params["step_value"] = 1
        print('step_value=',params["step_value"])
    elif event.key == '0':  # reset to 0th frame/file
        params["current_idx"] = 0
        plot_substrate()
    else:
        print('press', event.key)
def visualize_microenvironment(field_index):
    try:
      # apparently we need mpl's Qt backend to do keypresses
    #  matplotlib.use("Qt5Agg")
      matplotlib.use("TkAgg")
    except:
      print("\n---Error: cannot use matplotlib's TkAgg backend")
    #  print("Consider installing Anaconda's Python 3 distribution.")
      raise
#    current_idx = 0
#    frame_idx = 0
    #import matplotlib.pyplot as plt
#    if (len(sys.argv) < 3):
#      frame_idx = 0
#      field_index = 4
#    else:
#      kdx = 1
#      frame_idx = int(sys.argv[kdx])
#      kdx += 1
#      field_index = int(sys.argv[kdx])
    fig = plt.figure(figsize=(7,5.8))
    #ax = fig.gca()
    global params
    params = {}
    params["time_delay"] = 0.1
    params["count"] = -1
    params["cbar"] = None
    params["step_value"] =1
    params["current_idx"] = 0
    params["frame_idx"] = 0
    params["field_index"] = field_index
#    global time_delay
#    time_delay = 0.1
#    count = -1
#
#    cbar = None
#
#    #def plot_substrate(FileId):
#
#    step_value = 1

#    params = plot_substrate(params)
    plot_substrate()
    print("\nNOTE: click in plot window to give it focus before using keys.")
    fig.canvas.mpl_connect('key_press_event', press)
    #plot_substrate(frame_idx)
    plt.show()