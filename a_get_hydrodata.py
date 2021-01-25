# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:30:42 2020

Functions to set up hydro raster files for use in Python.

@author: clhenry
"""

import os
import numpy as np
from import_raw_tiffs import import_raw_tiff


# Save as numpy arrays of catchment area, slope, flow:
def get_hydrodata(folder, catcharea_infile, slope_infile, flow_infile, catcharea_outfile, slope_outfile, flow_outfile):
    
    wd = os.getcwd()
    folder_hydro = folder + '/Hydro'
    
    # Get raw data as numpy array
    catchment_area = import_raw_tiff(folder_hydro, catcharea_infile)
    slope = import_raw_tiff(folder_hydro, slope_infile)
    mean_flow = import_raw_tiff(folder_hydro, flow_infile)
    
    # Save files
    area_file = wd + '/Outputs/' + folder + '/' + catcharea_outfile
    np.save(area_file, catchment_area)

    slope_file = wd + '/Outputs/' + folder + '/' + slope_outfile
    np.save(slope_file, slope)

    flow_file = wd + '/Outputs/' + folder + '/' + flow_outfile
    np.save(flow_file, mean_flow)
