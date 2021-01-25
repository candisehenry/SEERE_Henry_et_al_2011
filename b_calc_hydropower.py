# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17

@author: hpaliwal, edited by candisehenry
"""

import numpy as np

# Hydropower density
def calc_hydropower(catchment_area, slope, mean_flow, step):

    '''
    Calculate the power density for locations with:
        catchment area (catchment_area) > min_catchment_area 100 km^2
        river gradient (slope) > min_slope 1%
        mean annual stream flow (mean_flow) > mean annual flow 4 m^3/s
    Output:
        power_density = hydro power density in kW/km^2
    '''
    
### Constants:

    # Hydraulic efficiency of turbine
    hydraulic_eff = 0.9

    # Catchment size (used as minimum sizing criteria)
    #min_catchment_area = 0.1 # 10 # km^2

    # Mean annual flow rate (used as minimum sizing criteria) (https://www.sciencedirect.com/science/article/pii/S1364032114003967)
    min_mean_flow = 2 # m^3/s # must be same as in costhydro()

    # River gradient percent (used as minimum sizing criteria) (https://www.sciencedirect.com/science/article/pii/S1364032114003967)
    min_slope = 0.5 / 100 # 0.5%
    
    # Upper bound of power output
    max_power_density = 10**9 # kW (=1 GW)
    
    # Density of water
    rho = 1000 # kg/m^3

    # Acceleration due to gravity
    g = 9.81 # m/s^2

    # Get raster info
    x = step
    cell_width = x * 111 * 1000 # m # 111 km is ~distance of 1 degree at 15 degrees latitude

### Power Density:

    # Select areas with slope greater than min_slope
    slope[slope < min_slope] = 0 

    # Convert slope to height difference
    height = np.multiply(slope, cell_width) # m
    
    # Eliminate flows that don't meet Szabo criteria
    mean_flow[mean_flow < min_mean_flow] = 0

    # Calculate power density per catchment area (should be same for all cells in catchment)
    power_density = (hydraulic_eff * rho * g * np.multiply(mean_flow, height)) / 1000 / catchment_area # kW/km^2

    # Remove data errors
    power_density[np.isinf(power_density)] = 0

    # Set power_max where power is greater than upper bound
    power_density[power_density > max_power_density] = max_power_density
    
    # Remove data where insufficient catchment area
    #power_density[catchment_area < min_catchment_area] = 0
    
    return power_density
