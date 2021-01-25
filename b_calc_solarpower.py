# -*- coding: utf-8 -*-
"""
Created on Mon May 18 2020

@author: hpaliwal, edited by candisehenry
"""

import os
import numpy as np
from a_get_solarwinddata import import_raw_nasa_data

def calc_solar(folder, infile, outfile1, outfile2, run):

    # Get raw NASA data file
    wd = os.getcwd()
    input_file = wd + '/Data/' + folder + '/' + infile
    rawdat = import_raw_nasa_data(input_file)
    
    if run == 'mean':
        # Get temperature data (average across years)
        temp_dat = rawdat.loc[rawdat['PARAMETER'] == 'T2M']
        temp = temp_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).mean()
        # Get insolation on a horizontal surface data (average across years)
        ihs_dat = rawdat.loc[rawdat['PARAMETER'] == 'ALLSKY_SFC_SW_DWN'] # kWh/m^2/day
        ihs = ihs_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).mean()
    elif run == 'max':
        temp_dat = rawdat.loc[rawdat['PARAMETER'] == 'T2M_MAX']
        temp = temp_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).max()
        ihs_dat = rawdat.loc[rawdat['PARAMETER'] == 'ALLSKY_SFC_SW_DWN'] # kWh/m^2/day
        ihs = ihs_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).min()
    elif run == 'min':
        temp_dat = rawdat.loc[rawdat['PARAMETER'] == 'T2M_MIN']
        temp = temp_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).min()
        ihs_dat = rawdat.loc[rawdat['PARAMETER'] == 'ALLSKY_SFC_SW_DWN'] # kWh/m^2/day
        ihs = ihs_dat.groupby(['LAT', 'LON', 'PARAMETER'], as_index=False).max()
    
    # Set input parameters
    latmin = min(ihs['LAT'])
    latmax = max(ihs['LAT'])
    longmin = min(ihs['LON'])
    longmax = max(ihs['LON'])
    
    # Constants that should not be changed
    y_pv =  1 # y_pv is the rated capacity of the PV array, meaning its power output under standard test conditions [kW/1-kW array]
    gt_stc = 1 # gt_stc is the solar radiation incident on the PV array in the current time step at standard test conditions [1 kW/m2] 
    alpha_p = -0.48 # alpha_p is the temperature coefficient of power [%/°C] 
    temp_c_stc = 25 # temp_c_stc is the PV cell temperature in the current time step under standard test conditions [25 °C] 
    
    # Preallocate values and matrix
    a,_ = np.shape(ihs)
    pvpwr = []
    solarcf = []
    
    # Go through all lines in data file and calculate only for the requested polygon
    for x in range(a):
        if ihs.loc[x,'LAT'] >= latmin and ihs.loc[x,'LAT'] <= latmax and ihs.loc[x,'LON'] >= longmin and ihs.loc[x,'LON'] <= longmax:
            
            pvlist = []
            cflist = []
            
            # Add lat and lon
            pvlist.append(ihs.loc[x,'LAT'])
            pvlist.append(ihs.loc[x,'LON'])

            cflist.append(ihs.loc[x,'LAT'])
            cflist.append(ihs.loc[x,'LON'])
            
            for y in range(4, 4+13): # 12 months and annual average
                # Calculate PV power
                gt_avg = ihs.iloc[x,y] / 24 # kW/m^2
                temp_c = temp.iloc[x,y]
                pvpower = y_pv * gt_avg / gt_stc * (1 + alpha_p / 100.0 * (temp_c - temp_c_stc)) # kW/1-kW array
                
                # Calculate capacity factor
                cf = gt_avg / gt_stc
                
                # Add pv power corresponding to that month
                pvlist.append(pvpower)
                cflist.append(cf)
                
            if np.size(pvpwr) == 0:
                pvpwr = pvlist
                solarcf = cflist
            else:
                pvpwr = np.vstack((pvpwr,pvlist))
                solarcf = np.vstack((solarcf,cflist))
    
    # Save file
    pv_file = open(wd + '/Outputs/' + folder + '/' + outfile1, 'w')
    np.savetxt(pv_file, pvpwr, fmt='%8.4f')
    pv_file.close()
    
    cf_file = open(wd + '/Outputs/' + folder + '/' + outfile2, 'w')
    np.savetxt(cf_file, solarcf, fmt='%8.4f')
    cf_file.close()