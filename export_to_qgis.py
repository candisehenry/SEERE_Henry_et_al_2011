# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 17:06:02 2020

@author: clhenry
"""
'''
from osgeo import ogr, gdal, osr
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import * 
from shapely.geometry import shape
import rasterio.features
from matplotlib.patches import Polygon
import fiona
import matplotlib.pyplot as plt
'''

import os
import numpy as np

def npy_to_geotiff(folder, infile, outfile, coord, step, nrows, ncols):

    # Import numpy files of LCOE
    wd = os.getcwd()
    
    file = wd + '/Outputs/' + folder + '/' + infile
    
    numpy_file = np.load(file, allow_pickle=True)
    
    # Export files as 
    out = wd + '/Outputs/' + folder + '/' + outfile
    f = open(out, 'w')
    f.write("NCOLS " + str(ncols) + "\n")
    f.write("NROWS " + str(nrows) + "\n")
    f.write("XLLCORNER " + str(coord[2]) + "\n")
    f.write("YLLCORNER " + str(coord[0]) + "\n")
    f.write("CELLSIZE " + str(step) + "\n")
    f.write("NODATA_VALUE " + str(-9999) + "\n")
    np.savetxt(f, numpy_file)
    f.close()
    
    
    
'''
# Import Guatemala country shapefile
shp_path = wd + '/Data/Guatemala/GTM_adm/'
shape = fiona.open(shp_path + 'GTM_adm0.shp')

# Create polygon from shapefile
pol = shape.next()
geom = pol['geometry']
poly_data = pol["geometry"]["coordinates"][0]
poly = Polygon(poly_data)


# Set up meshgrid
ihsfile = wd + '/Outputs/guatemala_pv_power.txt'
ihs = np.genfromtxt(ihsfile)

lat_min = min(ihs[:,0])
lat_max = max(ihs[:,0])
long_min = min(ihs[:,1])
long_max = max(ihs[:,1])
step = 1.0 / 120.0
x = np.arange(long_max, long_min, -step)
y = np.arange(lat_min, lat_max, step)
grid_x, grid_y = np.meshgrid(x,y)

# Save figures as raster / .png
CS = plt.pcolor(grid_x, grid_y, solar)
plt.colorbar()
plt.savefig(wd + '/Outputs/wind_linear.png', dpi = 300)
plt.show()
plt.close()
'''
