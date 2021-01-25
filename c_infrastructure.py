# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:39:33 2020

Create existing infrastructure file.

@author: clhenry
"""

import os
import numpy as np
from scipy.interpolate import griddata
from import_raw_tiffs import import_raw_tiff

# Collect all infrastructure tiffs into one numpy array
def get_infrastructure(folder, infile, outfile, step_from, step_to, coord_from, coord_to):

    # infile should be list with (1) transmission and (2) urban areas filenames
    
    wd = os.getcwd()

    data_new = []

    # Get grid for correct coordinates range (to)
    lat_min_to = coord_to[0]
    lat_max_to = coord_to[1]
    long_min_to = coord_to[2]
    long_max_to = coord_to[3]
    x = np.arange(long_min_to, long_max_to, step_to)
    y = np.arange(lat_max_to, lat_min_to, -step_to)
    grid_x, grid_y = np.meshgrid(x,y)

    for i in range(len(infile)):
        
        fol = folder[i+1]
        inf = infile[i]
        step_f = step_from[i]
        coord_f = coord_from[i]
    
        data = import_raw_tiff(fol, inf)
        
        # Reshape dataset 1 with existing coordinates range (from)
        lat_min_from = coord_f[0]
        lat_max_from = coord_f[1]
        long_min_from = coord_f[2]
        long_max_from = coord_f[3]
        x_from = np.arange(long_min_from, long_max_from, step_f)
        y_from = np.arange(lat_max_from, lat_min_from, -step_f)
        grid_x_from, grid_y_from = np.meshgrid(x_from, y_from)
        
        lon_from = grid_x_from.reshape(-1,1, order='F')
        lat_from = grid_y_from.reshape(-1,1, order='F')
        data = data.reshape(-1,1, order='F')
        points = np.concatenate([lon_from, lat_from, data], axis=1)
    
        # Remove values from existing data (from) that are outside coordinate range of new data (from)
        points = points[points[:,0] > long_min_to]
        points = points[points[:,0] < long_max_to]
        points = points[points[:,1] > lat_min_to]
        points = points[points[:,1] < lat_max_to]  
    
        # Regrid
        data = griddata(points[:,0:2], points[:,2], (grid_x, grid_y), method = 'nearest')
        
        # Combine infrastructure and urban areas dataset into one
        data_new =+ data

    # Turn tiff into black and white
    data_new[data_new > 0] = 1
    
    # Save file as numpy array
    infr_file = wd + '/Outputs/' + folder[0] + '/' + outfile
    np.save(infr_file, data_new)


# Collect all infrastructure tiffs into one numpy array (less clean code)
def get_infrastructure2(folder, infile, outfile, step_from, step_to, coord_from, coord_to):

    # infile should be list with (1) transmission and (2) urban areas filenames
    
    wd = os.getcwd()

    folder1 = folder + '/Transmission/'
    infile1 = infile[0]
    step_from1 = step_from[0]
    coord_from1 = coord_from[0]

    folder2 = folder + '/LandCover/'
    infile2 = infile[1]
    step_from2 = step_from[1]
    coord_from2 = coord_from[1]
    
    data1 = import_raw_tiff(folder1, infile1)
    data2 = import_raw_tiff(folder2, infile2)
    
    # Reshape dataset 1 with existing coordinates range (from)
    lat_min_from1 = coord_from1[0]
    lat_max_from1 = coord_from1[1]
    long_min_from1 = coord_from1[2]
    long_max_from1 = coord_from1[3]
    x_from1 = np.arange(long_min_from1, long_max_from1, step_from1)
    y_from1 = np.arange(lat_max_from1, lat_min_from1, -step_from1)
    grid_x_from1, grid_y_from1 = np.meshgrid(x_from1, y_from1)
    
    lon_from1 = grid_x_from1.reshape(-1,1, order='F')
    lat_from1 = grid_y_from1.reshape(-1,1, order='F')
    data1 = data1.reshape(-1,1, order='F')
    points1 = np.concatenate([lon_from1, lat_from1, data1], axis=1)

    # Reshape dataset 2 with existing coordinates range (from)
    lat_min_from2 = coord_from2[0]
    lat_max_from2 = coord_from2[1]
    long_min_from2 = coord_from2[2]
    long_max_from2 = coord_from2[3]
    x_from2 = np.arange(long_min_from2, long_max_from2, step_from2)
    y_from2 = np.arange(lat_max_from2, lat_min_from2, -step_from2)
    grid_x_from2, grid_y_from2 = np.meshgrid(x_from2, y_from2)
    
    lon_from2 = grid_x_from2.reshape(-1,1, order='F')
    lat_from2 = grid_y_from2.reshape(-1,1, order='F')
    data2 = data2.reshape(-1,1, order='F')
    points2 = np.concatenate([lon_from2, lat_from2, data2], axis=1)

    # Get grid for correct coordinates range (to)
    lat_min_to = coord_to[0]
    lat_max_to = coord_to[1]
    long_min_to = coord_to[2]
    long_max_to = coord_to[3]
    x = np.arange(long_min_to, long_max_to, step_to)
    y = np.arange(lat_max_to, lat_min_to, -step_to)
    grid_x, grid_y = np.meshgrid(x,y)

    # Remove values from existing data (from) that are outside coordinate range of new data (from)
    points1 = points1[points1[:,0] > long_min_to]
    points1 = points1[points1[:,0] < long_max_to]
    points1 = points1[points1[:,1] < lat_max_to]  
    points1 = points1[points1[:,1] > lat_min_to]

    points2 = points2[points2[:,0] > long_min_to]
    points2 = points2[points2[:,0] < long_max_to]
    points2 = points2[points2[:,1] < lat_max_to]  
    points2 = points2[points2[:,1] > lat_min_to]

    # Regrid
    data1 = griddata(points1[:,0:2], points1[:,2], (grid_x, grid_y), method = 'nearest')
    data2 = griddata(points2[:,0:2], points2[:,2], (grid_x, grid_y), method = 'nearest')    
    
    # Combine infrastructure and urban areas dataset into one
    data = data1 + data2

    # Turn tiff into black and white
    data[data > 0] = 1
    
    # Save file as numpy array
    infr_file = wd + '/Outputs/' + folder + '/' + outfile
    np.save(infr_file, data)
