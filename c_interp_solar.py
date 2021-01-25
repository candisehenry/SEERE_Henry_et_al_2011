# -*- coding: utf-8 -*-
"""
Created on May 18 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import numpy as np
from scipy.interpolate import griddata

def interp_solar(folder, infile, outfile, step, coord, run):
    # Get insolation on a horizontal surface file
    wd = os.getcwd()
    ihsfile = wd + '/Outputs/' + folder + '/' + infile
    ihs = np.genfromtxt(ihsfile)
    
    # Set input parameters
    lat_min = coord[0]
    lat_max = coord[1]
    long_min = coord[2]
    long_max = coord[3]
    x = np.arange(long_min, long_max, step)
    y = np.arange(lat_max, lat_min, -step)
    grid_x, grid_y = np.meshgrid(x,y)
    
    # Get annual average ihs and temp input parameters
    a = len(ihs)
    annual_ihs = np.zeros([int(a)])
    for k in range(a):
        if run == 'mean':
            annual_ihs[k] = np.mean(ihs[k, 2:-1])
        elif run == 'max':
            annual_ihs[k] = np.min(ihs[k, 2:-1]) # inverse because 'max' = max temp = min pv power
        elif run == 'min':
            annual_ihs[k] = np.max(ihs[k, 2:-1])
    
    # Reshape data and preallocate matrices
    a1, b1 = np.shape(ihs[:, 0:2])
    points = np.zeros([a1, b1])
    points[:,0] = ihs[:,1]
    points[:,1] = ihs[:,0]
    
    # Interpolate original data into grids using griddata
    ihs_interp = griddata(points, annual_ihs, (grid_x, grid_y), method = 'linear')
    
    # Fill NAs at edges
    mask = np.isnan(ihs_interp)
    ihs_interp[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), ihs_interp[~mask])

    # Save file
    ihs_file = wd + '/Outputs/' + folder + '/' + outfile
    np.save(ihs_file, ihs_interp)

