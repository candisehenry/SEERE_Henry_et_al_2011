# -*- coding: utf-8 -*-
"""
Created on May 19 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import numpy as np
from scipy.interpolate import griddata

def interp_wind(folder, infile, outfile, step, coord):
    # Get wind power density file
    wd = os.getcwd()
    windfile = wd + '/Outputs/' + folder + '/' + infile
    wpd = np.genfromtxt(windfile)
    
    # Set input parameters
    lat_min = coord[0]
    lat_max = coord[1]
    long_min = coord[2]
    long_max = coord[3]
    x = np.arange(long_min, long_max, step)
    y = np.arange(lat_max, lat_min, -step)
    grid_x, grid_y = np.meshgrid(x,y)
    
    # Get minimum (of all monthly averages) wind input parameter
    a = len(wpd)
    annual_win = np.zeros([int(a)])
    for k in range(a):
    	annual_win[k] = np.min(wpd[k, 2:])
    
    # Reshape data and preallocate matrices
    a1, b1 = np.shape(wpd[:, 0:2])
    points = np.zeros([a1, b1])
    points[:,0] = wpd[:,1]
    points[:,1] = wpd[:,0]
    
    # Interpolate original data into grids using griddata
    wpd_interp = griddata(points, annual_win, (grid_x, grid_y), method = 'linear')
    
    # Fill NAs at edges
    mask = np.isnan(wpd_interp)
    wpd_interp[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), wpd_interp[~mask])
    
    # Save file
    wpd_file = wd + '/Outputs/' + folder + '/' + outfile
    np.save(wpd_file, wpd_interp)
