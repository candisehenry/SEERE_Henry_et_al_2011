# -*- coding: utf-8 -*-
"""
Created on May 18 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import numpy as np
from scipy.interpolate import griddata

def interp_enedem(folder, infile, outfile, step, coord):
    # Load population density data
    wd = os.getcwd()
    enedemfile = wd + '/Outputs/' + folder + '/' + infile
    enedem = np.genfromtxt(enedemfile)
    
   # Set input parameters
    lat_min = coord[0]
    lat_max = coord[1]
    long_min = coord[2]
    long_max = coord[3]
    x = np.arange(long_min, long_max, step)
    y = np.arange(lat_max, lat_min, -step)
    grid_x, grid_y = np.meshgrid(x,y)

    # Get energy demand data
    energydemand = enedem[:,2]
    
    # Reshape data and preallocate matrices
    a1, b1 = np.shape(enedem[:, 0:2])
    points = np.zeros([a1, b1])
    points[:,0] = enedem[:,1]
    points[:,1] = enedem[:,0]
    
    # Interpolate original data into higher resolution data using griddata
    enedem_interp = griddata(points, energydemand, (grid_x, grid_y), method = 'linear')
    #enedem_interp[enedem_interp < 0] = 0
    
    # Fill NAs at edges
    mask = np.isnan(enedem_interp)
    enedem_interp[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), enedem_interp[~mask])
    
    # Save file
    enedem_file = wd + '/Outputs/' + folder + '/' + outfile
    np.save(enedem_file, enedem_interp)
