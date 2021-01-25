# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 11:15:29 2020

Run static model in parallel.

@author: clhenry
"""

import numpy
import pp
from core_model import costwind, costsolarbattery, costtrans, lcoe, parent_sol, parent_win, parent_hyd, parent_fossil


def run_static_model(global_params, step, ongrid_or_offgrid,
                     endem, solarava, windava, solarcf, windcf, infr, trans,
                     flow, catcharea, slope):
    
    ## Run model
    # Set initial parameters
    print("Setting up parallel computing.")
    a1, b1 = numpy.shape(endem)
    start = int(0)
    end = int(a1)
    parts = int(20)
    incr = int((end - start) / parts)
    
    # Create jobserver
    job_server = pp.Server()
    ncpus = 56
    job_server.set_ncpus(ncpus)
    
    # Submit a job which will calculate LCOEs in spatial "slices"
    print("Starting solar and wind ", job_server.get_ncpus(), " workers.")
    for index in range(parts):
        
        print("index = ", index, " of ", parts)
        
        # parent - the function
        # (starti, endi) - tuple with arguments for part_sum
        # (sizebatteryarray,costwind,lcoe,) - tuple with functions on which function parent depends
        # (energydemandpart,solaravapart,windavapart) - tuple with module names which must be imported before parent execution
        
        # Collect sliced inputs for solar and wind
        starti = start + index * incr
        endi = min(start + (index + 1) * incr, end)
        endempart = endem[starti:endi, 0:b1]
        solaravapart = solarava[starti:endi, 0:b1]
        windavapart = windava[starti:endi, 0:b1]
        solarcfpart = solarcf[starti:endi, 0:b1]
        windcfpart = windcf[starti:endi, 0:b1]
        infrpart = infr[starti:endi, 0:b1]
        
        # Run job
        f1 = job_server.submit(parent_sol, 
                               args=(global_params, endempart, solaravapart, solarcfpart, infrpart, trans, step),
                               depfuncs=(costsolarbattery, costtrans, lcoe,),
                               modules=('numpy', 'pandas', 'math', 'scipy.interpolate'))
        
        f2 = job_server.submit(parent_win, 
                               args=(global_params, endempart, windavapart, windcfpart, infrpart, trans, step),
                               depfuncs=(costwind, costtrans, lcoe,),
                               modules=('numpy', 'pandas', 'math', 'scipy.interpolate'))

        f3 = job_server.submit(parent_fossil, 
                               args=(global_params, endempart, infrpart, trans, step, ongrid_or_offgrid),
                               depfuncs=(costtrans, lcoe,),
                               modules=('numpy', 'pandas', 'math', 'scipy.interpolate'))
        
        if index == 0:
            lcsolar, pv_price = f1()
            #lcwind, wind_price = f2()
            #lcfossil = f3()
        else:
            lcsolari, pv_pricei = f1()
            lcwindi, wind_pricei = f2()
            lcfossili = f3()
            lcsolar = numpy.vstack((lcsolar, lcsolari))
            pv_price = numpy.vstack((pv_price, pv_pricei))
            lcwind = numpy.vstack((lcwind, lcwindi))
            wind_price = numpy.vstack((wind_price, wind_pricei))
            lcfossil = numpy.vstack((lcfossil, lcfossili))
    
    
    ## Run hydro as entire spatial extent because small (and thus fast) enough
    print("Starting hydro.")
    lchydro, hydro_price, hydrocf = parent_hyd(global_params, endem, flow, catcharea, slope, infr, trans, step)
    
    return lcsolar, lcwind, lchydro, lcfossil

