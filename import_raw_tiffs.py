# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 18:30:18 2020

Import tiff and geotiff rasters from QGIS.

@author: clhenry
"""

import os
from PIL import Image
import numpy as np


# Read raw tiff outputs from QGIS
def import_raw_tiff(folder, infile):
    
    wd = os.getcwd()
    input_file = wd + '/Data/' + folder + '/' + infile
    
    im = Image.open(input_file)
    
    data_array = np.array(im)
    
    return data_array
