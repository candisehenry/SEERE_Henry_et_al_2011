# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 13:16:16 2020

@author: clhenry
"""

import os
import numpy
import pandas
import time
import ast
from d_setup_guatemala import d_setup
from parallel_static_run import run_static_model
from export_to_qgis import npy_to_geotiff


start_time = time.time()

wd = os.getcwd()


# User input for file suffixes
print("Which streamflow scenario do you want to run? (i.e., filename string before '_flow')")
print("(Will be suffix in run name.)")
flow_suffix = str(input('Enter streamflow run name (min, mean, max): '))

print("\n")    
print("\nWhich temperature scenario do you want to run?")
print("(Will be suffix in run name.)")
temp_suffix = str(input('Enter temperature run (min, mean, max): '))

print("\n")
wind_mult = str(input('Enter the wind speed multiplication factor (0.78, 1, 1.21): '))
wind_multiplier = float(wind_mult)
if wind_multiplier > 1:
    wind_suffix = "max"
elif wind_multiplier < 1:
    wind_suffix = "min"
elif wind_multiplier == 1:
    wind_suffix = "mean"

print("\n")    
print("\nOn-grid or off-grid demand?")
ongrid_or_offgrid = str(input('Enter on or off: '))


# User input for whether to run data set-up
print("\nDo you want to run data set-up?")
print("A blank or unknown input will default to yes.")

run_d_setup = str(input('Enter y/n: '))

step, coord_solwin, coord_hyd = d_setup(run_d_setup, flow_suffix, temp_suffix, wind_multiplier)

print("\nFinished running data set-up.\n")

## Create dictionary of csv input parameters
# Open csv file
paramDF = pandas.read_csv(wd + '/Data/SEERE Cost Inputs/input_params_GLEDS.csv')
global_params = dict(zip(paramDF.key, paramDF.value))

# Convert dictionary into floats and lists (not strings)
for key in global_params:
    try:
        global_params[key] = int(global_params[key])
    except ValueError:
        global_params[key] = ast.literal_eval(global_params[key])

for key in global_params:
    if key == list:
        global_params[key] = numpy.array(global_params[key])


## Call in numpy files
# Get demand (kWh/day), wind, and solar files
endem = numpy.load(wd + '/Outputs/Guatemala/guatemala_demand_interp_survey.npy', allow_pickle=True) # kWh per day
solarava = numpy.load(wd + '/Outputs/Guatemala/guatemala_pv_interp.npy', allow_pickle=True) # kW
windava = numpy.load(wd + '/Outputs/Guatemala/guatemala_wind_interp.npy', allow_pickle=True) # kW
solarcf = numpy.load(wd + '/Outputs/Guatemala/guatemala_scf_interp.npy', allow_pickle=True)
windcf = numpy.load(wd + '/Outputs/Guatemala/guatemala_wcf_interp.npy', allow_pickle=True)

# Get demand (kWh/day, size of hydro data) and hydro files
flow_rate = numpy.load(wd + '/Outputs/Guatemala/flow.npy', allow_pickle=True) # kWh
catcharea = numpy.load(wd + '/Outputs/Guatemala/catcharea.npy', allow_pickle=True)
slope = numpy.load(wd + '/Outputs/Guatemala/slope.npy', allow_pickle=True)

# Existing infrastructure
infr = numpy.load(wd + '/Outputs/Guatemala/guatemala_urbanareas.npy', allow_pickle=True)
trans = numpy.load(wd + '/Outputs/Guatemala/guatemala_transmission.npy', allow_pickle=True)


# Run model
lcsolar, lcwind, lchydro, lcfossil = run_static_model(global_params, step, ongrid_or_offgrid,
                            endem, solarava, windava, solarcf, windcf, infr, trans, 
                            flow_rate, catcharea, slope)

print("\nTime elapsed: ", time.time() - start_time, "s")

