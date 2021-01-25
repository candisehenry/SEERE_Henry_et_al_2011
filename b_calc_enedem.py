
# -*- coding: utf-8 -*-
"""
Created on May 18 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from import_raw_tiffs import import_raw_tiff

# Calculate energy demand using csv input
def calc_energydemand_csv(folder, infile, outfile, step = None, coord = None):
    # Load raw population data
    wd = os.getcwd()
    folder_pop = folder + '/PopDensity'
    input_file = wd + '/Data/' + folder_pop + '/' + infile

    pop = pd.read_csv(input_file, header = None)
    
    # Energy demand per person per day (608 kWh/person/yr with 20% reserve capacity from worlddata.info)
    enedem = 608.0 / 365.0 * 1.2  # kWh/person/day
    
    # Create new table with coordinates
    enedemand = np.array(pop)
    
    # Calculate energy demand per day per grid
    enedemand[:,2] = np.multiply(pop.loc[:,2], enedem) # kWh/day

    # Save file
    demand_file = open(wd + '/Outputs/' + folder + '/' + outfile, 'w')
    np.savetxt(demand_file, enedemand, fmt='%8.4f')
    demand_file.close()
    
    # Then need to run output through interp_enedem.py


# Calculate energy demand using geotiff input
# Don't need to run interpolate in this case
def calc_energydemand_tif(folder, infile, outfile, step_from, step_to, coord_from, coord_to):
    # Load raw population data
    wd = os.getcwd()
    folder_pop = folder + '/PopDensity'
    pop = import_raw_tiff(folder_pop, infile) # fraction or person
    
    # Clean data
    if pop.dtype != 'float':
        pop = pop.astype('float')
    pop[pop <= 0.0005] = 'NaN'
    
    # Energy demand per person per day (81 kWh/person/yr with 20% reserve capacity from worlddata.info)
    # If pop is fractional population distribution, then multiply by 15,989,387 (total population) to get total demand
    # enedem = 81.0 / 365.0 * 1.2 * 15989387 # kWh/person/day
    enedem = 1.5 * 1000 * 24 # kWh/day # where first number is capacity expansion goal
    
    # Calculate energy demand per day per grid
    enedemand = np.multiply(pop, enedem) # kWh/day

    # Reshape data with existing coordinates range (from)
    lat_min_from = coord_from[0]
    lat_max_from = coord_from[1]
    long_min_from = coord_from[2]
    long_max_from = coord_from[3]
    x_from = np.arange(long_min_from, long_max_from, step_from)
    y_from = np.arange(lat_max_from, lat_min_from, -step_from)
    grid_x_from, grid_y_from = np.meshgrid(x_from,y_from)
    
    lon_from = grid_x_from.reshape(-1,1, order='F')
    lat_from = grid_y_from.reshape(-1,1, order='F')
    enedemand = enedemand.reshape(-1,1, order='F')
    points = np.concatenate([lon_from, lat_from, enedemand], axis=1)

    # Get grid for correct coordinates range (to)
    lat_min_to = coord_to[0]
    lat_max_to = coord_to[1]
    long_min_to = coord_to[2]
    long_max_to = coord_to[3]
    x = np.arange(long_min_to, long_max_to, step_to)
    y = np.arange(lat_max_to, lat_min_to, -step_to)
    grid_x, grid_y = np.meshgrid(x,y)

    # Remove values from existing data (from) that are outside coordinate range of new data (from)
    points = points[points[:,0] > long_min_to]
    points = points[points[:,0] < long_max_to]
    points = points[points[:,1] < lat_max_to]  
    points = points[points[:,1] > lat_min_to]

    # Regrid
    enedem_interp = griddata(points[:,0:2], points[:,2], (grid_x, grid_y), method = 'nearest')

    # Save file
    enedem_file = wd + '/Outputs/' + folder + '/' + outfile
    np.save(enedem_file, enedem_interp)
    
    # Can use output as is because already gridded.

