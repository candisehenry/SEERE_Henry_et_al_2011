# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 08:34:14 2020

Get min/max coordinates from input files.

@author: clhenry
"""

import os
import numpy as np
from a_get_solarwinddata import import_raw_nasa_data

# Check coordinates from raw NASA data file.
def set_coordinates(folder, raw_dat_infile):

    wd = os.getcwd()
    
    # Load solar and wind data for collecting lat and lon min/max
    infile = wd + '/Data/' + folder + '/' + raw_dat_infile
    rawdat = import_raw_nasa_data(infile)
    
    # Set input parameters
    lat_min = min(rawdat['LAT'])
    lat_max = max(rawdat['LAT'])
    long_min = min(rawdat['LON'])
    long_max = max(rawdat['LON'])
        
    return [lat_min, lat_max, long_min, long_max]
    

# Check coordinates from already cleaned pv and wind txt files.
def set_coordinates2(solar_infile, wind_infile):

    wd = os.getcwd()
    
    # Load solar and wind data for collecting lat and lon min/max
    ihsfile = wd + '/Outputs/' + solar_infile
    ihs = np.genfromtxt(ihsfile)
    
    windfile = wd + '/Outputs/' + wind_infile
    wpd = np.genfromtxt(windfile)

    # Set input parameters
    lat_min1 = min(ihs[:,0])
    lat_max1 = max(ihs[:,0])
    long_min1 = min(ihs[:,1])
    long_max1 = max(ihs[:,1])
    
    # Set input parameters
    lat_min2 = min(wpd[:,0])
    lat_max2 = max(wpd[:,0])
    long_min2 = min(wpd[:,1])
    long_max2 = max(wpd[:,1])
    
    if lat_min1 == lat_min2 and lat_max1 == lat_max2 and long_min1 == long_min2 and long_max1 == long_max2:
        
        return [lat_min1, lat_max1, long_min1, long_max1]
    
    else:
        
        print("Coordinates in solar and wind files don't match.")