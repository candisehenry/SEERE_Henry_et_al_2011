# -*- coding: utf-8 -*-
"""
Created on Mon Jun 1 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import numpy as np
from scipy.interpolate import griddata
from a_get_solarwinddata import import_raw_nasa_data
from import_raw_tiffs import import_raw_tiff

# Calculate wind power using csv input
def calc_wind_csv(folder, infile, outfile1, outfile2):

    # Get raw NASA data file
    wd = os.getcwd()
    input_file = wd + '/Data/' + folder + '/' + infile
    rawdat = import_raw_nasa_data(input_file)
    
    # Get insolation on a horizontal surface data (average across years)
    wind_dat = rawdat.loc[rawdat['PARAMETER'] == 'WS50M'] # m/s
    wind = wind_dat.groupby(['LAT', 'LON'], as_index=False).mean()
    
    # Get temperature data (average across years)
    temp_dat = rawdat.loc[rawdat['PARAMETER'] == 'T2M'] # C
    temp = temp_dat.groupby(['LAT', 'LON'], as_index=False).mean()
    
    # Set input parameters
    latmin = min(wind['LAT'])
    latmax = max(wind['LAT'])
    longmin = min(wind['LON'])
    longmax = max(wind['LON'])
    area = 1 # m^2
    rated_speed = 12 # m/s
    
    # Preallocate variables
    a, _ = np.shape(wind)
    windpwr = []
    windcf = []
    
    # Go through all lines in data file and calculate only for the requested polygon
    for x in range(a):
        if wind.loc[x,'LAT'] >= latmin and wind.loc[x,'LAT'] <= latmax and wind.loc[x,'LON'] >= longmin and wind.loc[x,'LON'] <= longmax:
    
            windlist = []
            cflist = []
    
            # Add lat and lon
            windlist.append(wind.loc[x,'LAT'])
            windlist.append(wind.loc[x,'LON'])

            cflist.append(wind.loc[x,'LAT'])
            cflist.append(wind.loc[x,'LON'])
    
            for y in range(3, 3+13): # 12 months and annual average
                # Calculate wind power
                wind_speed = wind.iloc[x,y]
                air_temp = temp.iloc[x,y]

                # Convert air temperature into average air density
                air_density = -0.003780 * float(air_temp) + 1.286400 # kg/m^3

                # Calculate wind power capacity for turbine with swept area = 1 m^2
                if wind_speed < 3:
                    windpower = 0
                elif wind_speed >= 3 and wind_speed < 24:
                    windpower = 0.5 * area * (wind_speed**3) * air_density * 1/1000 # kW (1 kg*m^2/s^3 = 1 watt)
                else:
                    windpower = 0

                # Calculate capacity factor (https://pubs.rsc.org/en/content/articlelanding/2018/ee/c7ee03029k#!divAbstract)
                if wind_speed < 3:
                    cf = 0
                elif wind_speed >= 3 and wind_speed < 12:
                    cf = wind_speed**3 / rated_speed**3
                elif wind_speed >= 12 and wind_speed < 24:
                    cf = 1
                else:
                    cf = 0

                # Add wind power corresponding to that month
                windlist.append(windpower)
                cflist.append(cf)
                
            if np.size(windpwr) == 0:
                windpwr = windlist
                windcf = cflist
            else:
                windpwr = np.vstack((windpwr,windlist))
                windcf = np.vstack((windcf,cflist))
    
    # Save file
    wind_file = open(wd + '/Outputs/' + folder + '/' + outfile1, 'w')
    np.savetxt(wind_file, windpwr, fmt='%8.4f')
    wind_file.close()
    
    cf_file = open(wd + '/Outputs/' + folder + '/' + outfile2, 'w')
    np.savetxt(cf_file, windcf, fmt='%8.4f')
    cf_file.close()
    
    # Then need to run output through interp_wind.py


# Calculate wind power using geotiff input
# Don't need to run interpolate in this case
def calc_wind_tif(folder, infile_wndspd, infile_pwrdns, outfile1, outfile2, 
                  step_from, step_to, coord_from, coord_to, variability):
    
    # Get tiffs
    wd = os.getcwd()
    folderpath = folder + '/Wind/'
    power_density = import_raw_tiff(folderpath, infile_pwrdns)
    wind_speed = import_raw_tiff(folderpath, infile_wndspd)

    # Set constants
    rated_speed = 12 # m/s

    # Preallocate matrices and get data size
    rows, cols = wind_speed.shape
    cf = np.zeros([rows, cols])

    # Convert power density from W/m^2 to kW (for turbine with swept area = 1 m^2)
    power_density[power_density < 1] = 0
    power_density = power_density * 1/1000
    
    # Calculate capacity factor (https://pubs.rsc.org/en/content/articlelanding/2018/ee/c7ee03029k#!divAbstract)
    for i in range(rows):
        for j in range(cols):
            # Set wind variability
            wind_speed[i,j] = wind_speed[i,j] * variability

            if wind_speed[i,j] < 3:
                cf[i,j] = 0
            elif wind_speed[i,j] >= 3 and wind_speed[i,j] < 12:
                cf[i,j] = wind_speed[i,j]**3 / rated_speed**3
            elif wind_speed[i,j] >= 12 and wind_speed[i,j] < 24:
                cf[i,j] = 1
            else:
                cf[i,j] = 0
    
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
    power_density = power_density.reshape(-1,1, order='F')
    cf = cf.reshape(-1,1, order='F')
    points1 = np.concatenate([lon_from, lat_from, power_density], axis=1)
    points2 = np.concatenate([lon_from, lat_from, cf], axis=1)

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
    power_density_interp = griddata(points1[:,0:2], points1[:,2], (grid_x, grid_y), method = 'nearest')
    cf_interp = griddata(points2[:,0:2], points2[:,2], (grid_x, grid_y), method = 'nearest')

    # Save files
    wpd_file = wd + '/Outputs/' + folder + '/' + outfile1
    np.save(wpd_file, power_density_interp)

    wcf_file = wd + '/Outputs/' + folder + '/' + outfile2
    np.save(wcf_file, cf_interp)
    
    # Can use output as is because already gridded.


# Calculate wind power using wind speed geotiff input and temperature npy input
# Don't need to run interpolate in this case
def calc_wind_tif_v2(folder, infile_wndspd, infile_temp, outfile1, outfile2, 
                  step_from_wndspd, step_from_temp, step_to, 
                  coord_from_wndspd, coord_from_temp, coord_to, 
                  variability):
    
    # Get tiffs
    wd = os.getcwd()
    folderpath = folder + '/Wind/'
    wind_speed = import_raw_tiff(folderpath, infile_wndspd)

    # Get temperature numpy
    air_temp = np.load(folderpath + infile_temp, allow_pickle=True) # kWh
    
    # Set constants
    area = 1 # m^2
    rated_speed = 12 # m/s
    
    # Get grids to set up regridding
    lat_min_from = coord_from_wndspd[0]
    lat_max_from = coord_from_wndspd[1]
    long_min_from = coord_from_wndspd[2]
    long_max_from = coord_from_wndspd[3]
    x_from = np.arange(long_min_from, long_max_from, step_from_wndspd)
    y_from = np.arange(lat_max_from, lat_min_from, -step_from_wndspd)
    grid_x_from_wndspd, grid_y_from_wndspd = np.meshgrid(x_from,y_from)
    
    lat_min_from = coord_from_temp[0]
    lat_max_from = coord_from_temp[1]
    long_min_from = coord_from_temp[2]
    long_max_from = coord_from_temp[3]
    x_from = np.arange(long_min_from, long_max_from, step_from_temp)
    y_from = np.arange(lat_max_from, lat_min_from, -step_from_temp)
    grid_x_from_temp, grid_y_from_temp = np.meshgrid(x_from,y_from)
    
    # Reshape data with existing coordinates range (from)
    lon_from = grid_x_from_wndspd.reshape(-1,1, order='F')
    lat_from = grid_y_from_wndspd.reshape(-1,1, order='F')
    wind_speed = wind_speed.reshape(-1,1, order='F')
    points1 = np.concatenate([lon_from, lat_from, wind_speed], axis=1)

    lon_from = grid_x_from_temp.reshape(-1,1, order='F')
    lat_from = grid_y_from_temp.reshape(-1,1, order='F')
    air_temp = air_temp.reshape(-1,1, order='F')
    points2 = np.concatenate([lon_from, lat_from, air_temp], axis=1)

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
    wind_speed_regrid = griddata(points1[:,0:2], points1[:,2], (grid_x, grid_y), method = 'linear')
    air_temp_regrid = griddata(points2[:,0:2], points2[:,2], (grid_x, grid_y), method = 'linear')

    # Preallocate matrices and get data size
    rows, cols = wind_speed_regrid.shape
    cf = np.zeros([rows, cols])
    power_density = np.zeros([rows, cols])

    # Calculate CF and power density
    for i in range(rows):
        for j in range(cols):
            # Set wind variability
            wind_speed_regrid[i,j] = wind_speed_regrid[i,j] * variability

            # Convert air temperature into average air density
            air_density = -0.003780 * float(air_temp_regrid[i,j]) + 1.286400 # kg/m^3

            # Calculate wind power capacity for turbine with swept area = 1 m^2
            if wind_speed_regrid[i,j] < 3:
                power_density[i,j] = 0
            elif wind_speed_regrid[i,j] >= 3 and wind_speed_regrid[i,j] < 24:
                power_density[i,j] = 0.5 * area * (wind_speed_regrid[i,j]**3) * air_density * 1/1000 # kW (1 kg*m^2/s^3 = 1 watt)
            else:
                power_density[i,j] = 0

            # Calculate capacity factor (https://pubs.rsc.org/en/content/articlelanding/2018/ee/c7ee03029k#!divAbstract)
            if wind_speed_regrid[i,j] < 3:
                cf[i,j] = 0
            elif wind_speed_regrid[i,j] >= 3 and wind_speed_regrid[i,j] < 12:
                cf[i,j] = wind_speed_regrid[i,j]**3 / rated_speed**3
            elif wind_speed_regrid[i,j] >= 12 and wind_speed_regrid[i,j] < 24:
                cf[i,j] = 1
            else:
                cf[i,j] = 0

    # Save files
    wpd_file = wd + '/Outputs/' + folder + '/' + outfile1
    np.save(wpd_file, power_density)

    wcf_file = wd + '/Outputs/' + folder + '/' + outfile2
    np.save(wcf_file, cf)
    
    # Can use output as is because already gridded.
