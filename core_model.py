# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:03:48 2020

Cleaned up version of pp_par_sol1_3d.py by hpaliwal

Equations to be called by run_SEACART.py

@author: clhenry
"""

import numpy
import math
import scipy.interpolate

# Total (sized) wind cost based on power density data ($/kW)
def costwind(energy_demand, power_avail, cf):
    
    # Input parameters
    eff = 0.3   # efficiency of wind to electricity; Betz's law of max power coefficient (https://www.raeng.org.uk/publications/other/23-wind-turbine)
    blades = 3  # number of blades per turbine

    # Peak wind hours per day
    peak_wind_hrs = cf * 24

    # Blade radius
    needed_sizing = (energy_demand / (power_avail * peak_wind_hrs * eff * math.pi))**(1/2) # meters

    if needed_sizing > 70: # maximum turbine radius
        sizing = 70
        multiples = math.ceil(needed_sizing / 70) # number of turbines needed, after reaching max build more turbines not bigger turbines
    elif needed_sizing <= 15:
        sizing = 15
        multiples = 1
    else:
        sizing = needed_sizing
        multiples = 1

    # Get the multiples of needed turbine size to actual turbine size
    mult_sizing = multiples * sizing / needed_sizing

    # $ per installation based on diameter of blades (empirically-based costs from https://www.nrel.gov/docs/fy07osti/40566.pdf)
    blade_frac_of_turbine_cost = 0.197 # (source: https://www.nrel.gov/docs/fy18osti/72167.pdf)
    turbine_frac_of_total_cost = 0.679
    pricing = multiples * blades * (2.5215203120 * (sizing)**2.7865867031) / blade_frac_of_turbine_cost / turbine_frac_of_total_cost # 2.28 * (0.2106 * sizing**2.6578) from NREL doc above
    #pricing = multiples * blades * (2.5215203120 * (sizing)**2.7865867031) / (0.2 * 0.64) # numbers from original SEACART

    if numpy.isfinite(pricing) == False:
        pricing = 0

    # # Price adjustment parameter
    # adjpricing = 1.0187598973 * (power_demand**(-0.0386783672)) # dimensionless
    # 
    # if adjpricing > 1.163:
    #     adjpricing = 1.163
    # 
    # if adjpricing < 0.85:
    #     adjpricing = 0.85

    # Price per kW of wind
    total_rated_installed_power = mult_sizing * energy_demand / peak_wind_hrs # kW # equals power_demand/CF
    windpriceperkw = pricing / total_rated_installed_power # $/kW
    # windpriceperkw = adjpricing * pricing / power_demand # $/kW # adjpricing from original SEACART

    if numpy.isfinite(windpriceperkw) == False:
        windpriceperkw = 0 
    
    return windpriceperkw


# Solar cost interpolation (for entire solar system, not just PV modules)
def interpscale(x):

    '''
    Input PV size as x (kW) to output interpolated price ($/kW)
    '''

    # Solar panel system size kW
    pv_size_in = numpy.array([0.1,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,30.,50.,100.,90000000.])
    
    # Price per $/kW
    price_in = 1 / 3.04 * numpy.array([5.0,4.47,4.25,4.05,3.88,3.73,3.61,3.50,3.41,3.33,3.27,3.22,3.19,3.15,3.13,3.11,3.09,3.08,3.06,3.04,3.04,3.04,3.04,3.04])
    
    # Function for price wrt scale
    pr_wrt_scl_f = scipy.interpolate.interp1d(pv_size_in, price_in)
    
    return pr_wrt_scl_f(x)


# Total (sized) solar and battery cost based on power density data ($/kW)
def costsolarbattery(energy_demand, power_avail, cf):

    '''
    Inputs:
        energy_demand = energy demand in kWh/day
        power_avail = solar availability in kW
        peak_sun_hrs = number of peak sun hours per day
        
    Outputs:
        bat_price = capital cost per kWh battery
        pv_price = capital cost per kW PV
        # capex = capital cost PV and battery combined
        # pv_capex_perkwh = capital cost per kWh PV
        # total_power_available = nameplate capacity of system
    '''
       
    # Peak sun hours per day
    peak_sun_hrs = cf * 24
    
##### Inputs that can be changed depending on the  battery and solar panel/module type:
    ## See: Trina Solar, Longi, Risen
    
    # Bus voltage for battery connection
    bbv = 24 # volts
    
    # Days of storage desired/required
    day_auto = 3 # days
    
    # Allowable depth of discharge limit
    dod = 0.8 # fraction
    
    # Price of battery: Crown CR220, 6V Flooded
    bat_price = 145 # $
    
    # Capacity of selected battery from battery specification sheet
    bat_cap_amp_hr_spec = 220 # amp-hours
    
    # Voltage of selected battery from battery specification sheet
    bat_volt_spec = 6 # volts
    
    # Battery roundtrip efficeincy
    bat_round_trp_eff = 0.85 # fraction
    
    # PV panel price (Trina Solar 72-CELL TALLMAX M DE14A(II), specs at https://static.trinasolar.com/sites/default/files/Datasheet_Tallmax_1500V%20M%20Plus_2019_May.pdf)
    pv_unit_price = 230 # $/panel
    
    # PV panel max power output at STC (Trina Solar 72-CELL TALLMAX M DE14A(II))
    pv_power_max = 380 # watts
    
    # PV panel derating factor
    derating_factor = 0.98

    # Contribution of the PV module to the total cost fraction (source: https://www.nrel.gov/docs/fy19osti/72133.pdf)
    PV_module_contribution = 0.5
        
##### Size battery array:
    
    # Total demand per day in amp-hours (where watt-hour = amp-hour * volt)
    amp_hr_Pday = (energy_demand * 1000) / bbv # amp-hours
    
    # Required battery charge capacity
    bat_cap_amp_hr = amp_hr_Pday * day_auto / dod # amp-hours
    
    # Number of batteries in parallel
    num_bat_par = math.ceil(bat_cap_amp_hr / bat_cap_amp_hr_spec)
    
    # Number of batteries in series
    num_bat_ser = math.ceil(bbv / bat_volt_spec)
     
    # Total number of batteries
    bat_tot = num_bat_par * num_bat_ser
    
##### Size PV array:

    # PV panel power output at STC (including battery roundtrip efficiency)
    pv_power_out = pv_power_max * derating_factor * bat_round_trp_eff * 1/1000 # kW
    
    # Total energy avail. for 1 panel of selected panel type = solar capacity of 1-kW panel * watt rating of selected panel type * number of hours
    tot_solar_ene_avail = power_avail * pv_power_out * peak_sun_hrs # kWh/day

    # Number of PV modules required to meet energy requirement
    num_pv_tot = math.ceil(energy_demand / tot_solar_ene_avail) # unitless
    # num_pv_tot = math.ceil( (energy_demand / bat_round_trp_eff) / tot_solar_ene_avail ) # unitless
    
    # Total PV installed power
    total_rated_installed_power = num_pv_tot * pv_power_out # kW
    
##### Price of battery and PV array:

    # Capex of total PV and battery setup
    # total_capex = num_pv_tot * pv_unit_price / PV_module_contribution + bat_tot * bat_price
    # pv_capex_perkWh =  num_pv_tot * pv_unit_price / PV_module_contribution / energy_demand # $/kWh
    
    # Battery price per kWh
    bat_price = bat_tot * bat_price / (total_rated_installed_power * 24) # $/kWh

    # Interpolated price of PV installation by converting total rated power
    #pv_price = interpscale(total_rated_installed_power) * (num_pv_tot * pv_unit_price / PV_module_contribution) / total_power_available # $/kW
    pv_price = (num_pv_tot * pv_unit_price / PV_module_contribution) / total_rated_installed_power # $/kW
    
    return bat_price, pv_price


# Total (sized) hydropower cost based on power density ($/kWh)
def costhydro(power_demand, slope, flow, power_avail, catch_area, step):

    '''
    Calculates the cost of installing a hydropower plant with a given power density and demand and net water head
    Uses correlations from and converts euros to dollars based on input
    
    Inputs:
        power_demand                        Power demand array (without raster info) (kW)
        slope                               Slope raster (in raster format without info) (%)
        flow                                Flow raster (in raster format without info) (m^3/s)
        power_avail                         Power density/available raster (in raster format without info) (kW/km^2)
        catch_area                          Area considered by the power_avail raster (km^2)
    
    Outputs:
        percentDemandMet
            percentDemandMetPelton          Percent that power density from hydropower can meet required demand by pelton turbines (in raster format without info)
            percentDemandMetFrancis         Percent that power density from hydropower can meet required demand by francis turbines (in raster format without info)
            percentDemandMetPropeller      Percent that power density from hydropower can meet required demand by kaplan and semi-kaplan turbines (in raster format without info)
            percentDemandMetKaplan          Percent that power density from hydropower can meet required demand by kaplan turbines (in raster format without info)
        Cmin                                Raster (in raster format without info) of lowest cost turbine price ($/kWh)
        turbineCosts
            pelton                          Raster (in raster format without info) of cost of pelton turbine ($/kWh)
            francis                         Raster (in raster format without info) of cost of francis turbine ($/kWh)
            propeller                       Raster (in raster format without info) of cost of propeller turbine ($/kWh)
            kaplan                          Raster (in raster format without info) of cost of kaplan turbine ($/kWh)
        turbinePower
            pelton                          Raster (in raster format without info) of energy output of pelton turbine (kWh)
            francis                         Raster (in raster format without info) of energy output of francis turbine (kWh)
            propeller                       Raster (in raster format without info) of energy output of propeller turbine (kWh)
            kaplan                          Raster (in raster format without info) of energy output of kaplan turbine (kWh)

    Costing from:
        Cost determination of the electro-mechanical equipment of a small hydro-power plant
        B. Ogayar*, P.G. Vidal, 2009
    '''

### Inputs:

    # Eliminate flows that don't meet Szabo criteria (https://www.sciencedirect.com/science/article/pii/S1364032114003967)
    mean_flow = 2 # m^3/s # minimum flowate to consider for streams (m^3/s), must be same as in hydroPower fun()
    if flow < mean_flow:
        flow = 0

    # Sizing and costing
    # (https://escholarship.org/content/qt0jb5v4df/qt0jb5v4df.pdf?t=oeis5t)
    # (https://energiatalgud.ee/img_auth.php/a/ab/Guide_on_How_to_Develop_a_Small_Hydropower_Plant.pdf)
    limitPowerOutputToDemand = True # Boolean: true limits the power output sizing to a maximum of the power demand, false has no maximum power output    
    standardPelton = [10000, 2, 1000, 0.9] # [[80000, 20, 1600, 0.9]] # array of rated standard conditions for pelton turbines  [[P, Q, H, efficiency], ...] (kW, m^3/s, m, fraction)
    standardFrancis = [10000, 20, 300, 0.9] # [[100000, 100, 300, 0.9]] # array of rated standard conditions for francis turbines  [[P, Q, H, efficiency], ...] (kW, m^3/s, m, fraction)
    standardKaplan = [10000, 50, 50, 0.9] # [[100000, 200, 70, 0.9]] # array of rated standard conditions for kaplan turbines  [[P, Q, H, efficiency], ...] (kW, m^3/s, m, fraction)
    #turbineFractionOfInstallationCost = 0.35 # fraction of installation costs accounted for by turbines

    # Inputs for interpolating turbine efficiencies (https://www.researchgate.net/figure/Efficiency-curve-of-various-turbine-designs-per-percentage-of-rated-flow-18_fig2_322639914)
    x = numpy.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    peltonEffCurve = numpy.array([0, 0.05, 0.6, 0.8, 0.89, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9])
    francisEffCurve = numpy.array([0, 0, 0, 0.3, 0.5, 0.6, 0.75, 0.85, 0.9, 0.9, 0.9])
    kaplanEffCurve = numpy.array([0, 0.4, 0.75, 0.85, 0.89, 0.9, 0.9, 0.89, 0.88, 0.87, 0.87])

    # Other inputs
    euro_to_dollar = 1.18 # $/euro

### Efficiency Interpolation:

    # Create interpolation functions    
    peltonInterpEff = scipy.interpolate.interp1d(x, peltonEffCurve)
    francisInterpEff = scipy.interpolate.interp1d(x, francisEffCurve)
    kaplanInterpEff = scipy.interpolate.interp1d(x, kaplanEffCurve)
    
    # Calculate percentage of standard flow for interpolation
    # Pelton
    percentStdFlowPelton = flow / standardPelton[1]
    if math.isnan(percentStdFlowPelton):
        percentStdFlowPelton = 0
    elif percentStdFlowPelton < 0:
        percentStdFlowPelton = 0
    else:
        percentStdFlowPelton = 1
    peltonEff = peltonInterpEff(percentStdFlowPelton)

    # Francis
    percentStdFlowFrancis = flow / standardFrancis[1]
    if math.isnan(percentStdFlowFrancis):
        percentStdFlowFrancis = 0
    elif percentStdFlowFrancis < 0:
        percentStdFlowFrancis = 0
    else:
        percentStdFlowFrancis = 1
    francisEff = francisInterpEff(percentStdFlowFrancis)

    # Kaplan
    percentStdFlowKaplan = flow / standardKaplan[1]
    if math.isnan(percentStdFlowKaplan):
        percentStdFlowKaplan = 0
    elif percentStdFlowKaplan < 0:
        percentStdFlowKaplan = 0
    else:
        percentStdFlowKaplan = 1
    kaplanEff = kaplanInterpEff(percentStdFlowKaplan)
    
### Sizing:

    # Prepare rasters
    cell_width = step * 111 * 1000 # m # 111km ~= 1 degree
    H = numpy.tan(numpy.deg2rad(slope * 100)) * cell_width # slope / 100 * cell_width
    if H <= 0:
        H = numpy.nan

    # Get power demand
    if power_demand <= 0 or power_avail <= 0:
        power_demand = numpy.nan

    # Get available power from power density
    power_avail = power_avail * catch_area
    if power_avail <= 0:
        power_avail = numpy.nan

    # Determine where demand is greater than energy available
    DemandGreaterThanEnergyAvailable = power_demand / power_avail
    if math.isnan(DemandGreaterThanEnergyAvailable):
        DemandGreaterThanEnergyAvailable = 0

### Pelton Size and Cost:
    
    ## Small hydro (<2 MW) cost data:
    # Cost data from Ogayar and Vidal 2009: https://www.sciencedirect.com/science/article/abs/pii/S0960148108001900
    # Open-access: https://dergipark.org.tr/tr/download/article-file/586584
    ## Larger hydro (>2 - 30 MW) cost data:
    # ORNL 2015: https://info.ornl.gov/sites/publications/files/Pub53978.pdf
    # ORNL 2015b: https://info.ornl.gov/sites/publications/files/Pub58666.pdf
    
    # Size power output relative to standard conditions for all turbines of each type [Pturbine = Ps*n*Q*H/(ns*Qs*Hs)]
    Ppelton = standardPelton[0] * (flow * H * peltonEff) / (standardPelton[1] * standardPelton[2] * standardPelton[3]) / 1000 # kW
    if Ppelton <= 0:
        Ppelton = numpy.nan

    if (limitPowerOutputToDemand == True):
        if Ppelton > power_demand:
            Ppelton = power_demand

    # Get number of turbines to build based on demand or power availability
    if DemandGreaterThanEnergyAvailable > 1:
        numPeltonTurbines = power_avail / Ppelton
    else:
        numPeltonTurbines = power_demand / Ppelton

    # Costing for power range below vs. above 2 MW
    if Ppelton < 2000:
        Cpelton = (17693 * Ppelton**(-0.3644725) * H**(-0.281735)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW
    else: 
        Cpelton = (17693 * 2000**(-0.3644725) * H**(-0.281735)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW 
        #Cpelton = (8717830 * (Ppelton/1000)**(0.975) * (H/3)**(-0.120)) / Ppelton # $/kW
    
    # Get percent of demand met by turbine generation
    percentDemandMetPelton = Ppelton / power_demand * 100
    
    # Replace bad values
    if math.isnan(numPeltonTurbines):
        numPeltonTurbines = 0
    elif math.isnan(percentDemandMetPelton):
        percentDemandMetPelton = 0
    elif math.isnan(Cpelton):
        Cpelton = 0
    elif math.isnan(Ppelton):
        Ppelton = 0

### Francis Size and Cost:

    # Size power output relative to standard conditions for all turbines of each type [Pturbine = Ps*n*Q*H/(ns*Qs*Hs)]
    Pfrancis = standardFrancis[0] * (flow * H * francisEff) / (standardFrancis[1] * standardFrancis[2] * standardFrancis[3]) / 1000 # kW
    if Pfrancis <= 0:
        Pfrancis = numpy.nan

    if (limitPowerOutputToDemand == True):
        if Pfrancis > power_demand:
            Pfrancis = power_demand

    # Get number of turbines to build based on demand or power availability
    if DemandGreaterThanEnergyAvailable > 1:
        numFrancisTurbines = power_avail / Pfrancis
    else:
        numFrancisTurbines = power_demand / Pfrancis

    # Costing for power range below vs. above 2 MW
    if Pfrancis < 2000:
        Cfrancis = (25698 * Pfrancis**(-0.560135) * H**(-0.127243)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW
    else:
        Cfrancis = (25698 * 2000**(-0.560135) * H**(-0.127243)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW
        #Cfrancis = (8717830 * (Pfrancis/1000)**(0.975) * (H/3)**(-0.120)) / Pfrancis # $/kW

    # Get percent of demand met by turbine generation
    percentDemandMetFrancis = Pfrancis / power_demand * 100
    
    # Replace bad values
    if math.isnan(numFrancisTurbines):
        numFrancisTurbines = 0
    elif math.isnan(percentDemandMetFrancis):
        percentDemandMetFrancis = 0
    elif math.isnan(Cfrancis):
        Cfrancis = 0
    elif math.isnan(Pfrancis):
        Pfrancis = 0
    
### Kaplan Size and Cost:

    # Size power output relative to standard conditions for all turbines of each type [Pturbine = Ps*n*Q*H/(ns*Qs*Hs)]
    PKaplan = standardKaplan[0] * (flow * H * kaplanEff) / (standardKaplan[1] * standardKaplan[2] * standardKaplan[3]) / 1000 # kW
    if PKaplan <= 0: 
        PKaplan = numpy.nan

    if (limitPowerOutputToDemand == True):
        if PKaplan > power_demand:
            PKaplan = power_demand

    # Get number of turbines to build based on demand or power availability
    if DemandGreaterThanEnergyAvailable > 1:
        numKaplanTurbines = power_avail / PKaplan
    else: 
        numKaplanTurbines = power_demand / PKaplan

    # Costing for power range below vs. above 2 MW
    if PKaplan < 2000:
        CKaplan = (19498 * PKaplan**(-0.58338) * H**(-0.113901)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW
    else:
        CKaplan = (19498 * 2000**(-0.58338) * H**(-0.113901)) * euro_to_dollar # / turbineFractionOfInstallationCost # $/kW
       # CKaplan = (8717830 * (PKaplan/1000)**(0.975) * (H/3)**(-0.120)) / PKaplan # $/kW

    # Get percent of demand met by turbine generation
    percentDemandMetKaplan = PKaplan / power_demand * 100

    # Replace bad values
    if math.isnan(numKaplanTurbines):
        numKaplanTurbines = 0
    elif math.isnan(percentDemandMetKaplan):
        percentDemandMetKaplan = 0
    elif math.isnan(CKaplan):
        CKaplan = 0
    elif math.isnan(PKaplan):
        PKaplan = 0

### Combine outputs to get least cost and percent demand met:

    #percentDemandMet = max(percentDemandMetPelton, percentDemandMetFrancis, percentDemandMetKaplan)

    if max(percentDemandMetPelton, percentDemandMetFrancis, percentDemandMetKaplan) == percentDemandMetPelton:
        Cmin = Cpelton             
        CF = percentStdFlowPelton
    
    if max(percentDemandMetPelton, percentDemandMetFrancis, percentDemandMetKaplan) == percentDemandMetFrancis:
        Cmin = Cfrancis
        CF = percentStdFlowFrancis
        
    if percentDemandMetFrancis == percentDemandMetPelton:
        Cmin = min(Cpelton, Cfrancis)
    
    if Cmin == Cpelton:
        CF = percentStdFlowPelton
    
    if max(percentDemandMetPelton, percentDemandMetFrancis, percentDemandMetKaplan) == percentDemandMetKaplan:
        Cmin = CKaplan
        CF = percentStdFlowKaplan
    
    if percentDemandMetKaplan == percentDemandMetPelton:
        Cmin = min(CKaplan, Cpelton)
        
    if Cmin == Cpelton:
        CF = percentStdFlowPelton
    
    if percentDemandMetKaplan == percentDemandMetFrancis:
        Cmin = min(CKaplan, Cfrancis)
    
    if Cmin == Cfrancis:
        CF = percentStdFlowFrancis

    # Convert NANs back to 0
    if math.isnan(CF):
        CF = 0
    elif CF > 1:
        CF = 1

    if math.isnan(Cmin):
        Cmin = 0
    elif Cmin > 10000:
        Cmin = 0

    return Cmin, CF


# Transmission cost ($)
def costtrans(i, j, trans_infr, step, capex):
    
    # Inputs:
    #   trans_infr          Numpy array of existing transmission infrastructure
    #   step                Cell size in degrees
    #   start_index         Get row index to ensure location of cell taken from correct 'slice'
    #   capex               Transmission cost in $/m
    # Outputs:
    #   transcost           Total transmission cost given transmission distance
    #   build_generation    New (reduced/sliced) energy demand matrix of all areas where generation will be built
    
    # Set constants
    km_len = 111.321 * step # km # length of 1 degree at the equator
    
    # Create masks
    # Mask where no demand but power available
    mask1 = (trans_infr > 0)
    row1, col1 = numpy.where(mask1 == True)

    # Get distances between current cell [i,j] and all cells with transmission
    array_dist = numpy.sqrt((i-row1)**2 + (j-col1)**2)

    # Find shortest distance to transmission line
    min_distance = min(array_dist)

    # Get transmission cost
    dist = km_len * min_distance
    transcost = dist * capex
    
    return transcost


# Calculate levelized cost of electricity
def lcoe(inputArray):

    # LCOE Calculations: http://large.stanford.edu/courses/2010/ph240/vasudev1/
    # Simplified LCOE Calculations: https://energy.utexas.edu/sites/default/files/UTAustin_FCe_LCOE_2016.pdf
    
    # Convert input array into variable names
    target_generation = inputArray[0]   # Target average production (kWh/day)
    cf = inputArray[1]                  # Average capacity factor (fraction of time at max)
    construction_time = inputArray[2]   # Total construction time - must match length of construction schedule array
    capital_cost = inputArray[3]        # Total nominal (overnight) busbar cost ($/kW)
    construction_sched = inputArray[4]  # Construction schedule array with fraction of total nominal (overnight) busbar cost spent over each year
    operation_years = inputArray[5]     # Number of years of operation
    operation_start = inputArray[6]     # Starting Year of operation (ex. 2020)
    depreciation_sched = inputArray[7]  # Depreciation schedule with percent depreciated each year (%)
    incr_capital_cost = inputArray[8]   # Incremental capital costs ($/kW/year)
    fixed_om_cost = inputArray[9]       # Fixed O&M costs ($/kW/year)
    var_om_cost = inputArray[10]        # Variable O&M costs ($/kWh)
    decommission_cost = inputArray[11]  # Decommisioning cost ($)
    inflation_rate = inputArray[12]     # Inflation rate (fraction)
    inflation_year = inputArray[13]     # Inflation base year (ex. 2016)
    storage_eff = inputArray[14]        # Storage efficiency (fraction)
    storage_cost = inputArray[15]       # Storage cost ($/kWh)
    storage_lifetime = inputArray[16]   # Storage lifetime in days (kWh)
    tax_rate = inputArray[17]           # Tax rate (%)
    debt_frac = inputArray[18]          # Debt fraction (fraction)
    debt_rate = inputArray[19]          # Debt rate (%)
    equity_rate = inputArray[20]        # Equity rate (%)
    interest_rate = inputArray[21]      # Interest rate (%) (tax/debt/equity rates and debt frac not needed if interest rate )
    transmission_cost = inputArray[22]  # Cost of transmission lines


    # Convert target output in kWh/day to kW
    target_output = target_generation / 24

    # Mask for construction or operational year (1 = operation, 0 = construction)
    operational_mask = numpy.append(numpy.zeros([construction_time]), numpy.ones([operation_years]))
    construction_mask = (operational_mask == 0) * 1

    # Create an array of years from 10 years before operation_start (construction) to operation_years after operation_start (operational)
    timeline = numpy.array(range(operation_start - construction_time, operation_start + operation_years))
    operation_period = timeline - operation_start + 1
    
    # Find nameplate capacity
    nameplate_capacity = target_output / cf # kW

    # Calculate inflation, weighted average cost of capital, and discount for each year
    inflation =  (1 + inflation_rate)**numpy.subtract(timeline, inflation_year,  dtype='float')
    if interest_rate == 0:
        wacc = ((1 - tax_rate / 100) * debt_frac * debt_rate) + ((1 - debt_frac) * equity_rate) # %
    else:
        wacc = interest_rate # %
    discount = numpy.power((1 + wacc / 100), operation_period) # (1+r)^t

    # Calculate construction costs
    construction_costs = numpy.append(construction_sched, numpy.zeros(operation_years)) * (capital_cost * nameplate_capacity + transmission_cost) * inflation # $

    # Calculate incremental capital costs and decommissioning
    incremental_cap_cost_and_decommiss = numpy.multiply(operational_mask, nameplate_capacity * incr_capital_cost) # $
    incremental_cap_cost_and_decommiss[len(operational_mask) - 1] += decommission_cost # $
    incremental_cap_cost_and_decommiss *= inflation

    # Calculate fixed operational and maintenance costs
    fixed_om_costs = operational_mask * nameplate_capacity * fixed_om_cost * inflation # $

    # Calculate consumed electricityType
    electricity_consumed = operational_mask * target_output * storage_eff * 8766 # kWh/year

    # Calculate variable O&M costs
    var_om_costs = var_om_cost * electricity_consumed * inflation # $/year

    # Battery replacements made every x years where x = years of operation / number of replacements rounded down
    index = numpy.where(timeline == operation_start - 1)[0][0]
    storage_replacements = int(operation_years * 365 / storage_lifetime) # multiply by 365 because battery lifetime in days
    replacement_gap = math.floor(operation_years / storage_replacements) # years between replacements
    replacement_idx = numpy.add(numpy.multiply(range(1, storage_replacements + 1), replacement_gap), index) # add index because only replace after operation starts

    # Calculate battery costs
    storage_costs = numpy.asarray((timeline == operation_start - 1), dtype='float')
    starting_storage_cap_inventory_cost = nameplate_capacity * 24 * storage_cost # $
    storage_costs[replacement_idx] = True
    storage_costs = numpy.multiply(storage_costs, starting_storage_cap_inventory_cost * inflation) # $

    # Total pre-depreciation costs
    total_predeprec_costs = numpy.add(construction_costs, numpy.add(incremental_cap_cost_and_decommiss, numpy.add(fixed_om_costs, numpy.add(var_om_costs, storage_costs)))) # $

    # Calculate discounted total pre-depreciation costs (NPV of pre-depreciation costs)
    discounted_total_predeprec_costs = total_predeprec_costs / discount # $ # NPV = initial investment / (1+r)^t

    # Calclulate inflated electricity consumed
    inflated_electricity_consumed = electricity_consumed * inflation # kWh/year

    # Calculate discounted inflated electricity consumed
    discounted_inflated_electricity_consumed = inflated_electricity_consumed / discount # kWh/year

    # Calculate pre-tax LCOE without depreciation tax shield
    pretax_LCOE = numpy.sum(discounted_total_predeprec_costs) / numpy.sum(discounted_inflated_electricity_consumed) # $/kWh

    # Calculate pre-tax taxable LCOE (pretax_LCOE - discounted construction costs)
    pretax_taxable_LCOE = pretax_LCOE - numpy.sum(construction_mask * discounted_total_predeprec_costs) / numpy.sum(discounted_inflated_electricity_consumed) # $/kWh

    # Calculate depreciation tax shield
    depreciation = numpy.sum(construction_costs) * operational_mask * numpy.append(numpy.zeros([construction_time]), numpy.array(depreciation_sched[:operation_years]) / 100) * tax_rate / 100 # $

    # Calculate discounted depreciation tax shield
    discounted_deprec = depreciation / discount # $

    # Calculate post-depreciation cost
    if (abs(tax_rate) > 0):
        postdeprec_cost = (numpy.sum(discounted_total_predeprec_costs) - numpy.sum(discounted_deprec) / (tax_rate / 100)) / numpy.sum(discounted_inflated_electricity_consumed) # $/kWh
        # postdeprec_tax_shield_LCOE = (numpy.sum(discounted_total_predeprec_costs) - numpy.sum(discounted_deprec)) / numpy.sum(discounted_inflated_electricity_consumed) # $/kWh

    else:
        postdeprec_cost = pretax_LCOE # $/kWh

    # Pretax profit
    pretax_profit = postdeprec_cost - pretax_taxable_LCOE # $/kWh
    
    # Calculate LCOE
    LCOE = pretax_LCOE + pretax_profit * (tax_rate / 100) / (1 - (tax_rate / 100)) # $/kWh

    return LCOE


# Use above functions to calculate LCOEs of solar
def parent_sol(p, ed, sa, scf, infr, trans, step):
    
    # Preallocate output matrices
    ae, be = numpy.shape(ed)
    lcs = numpy.zeros([ae, be])
    pv_price = numpy.zeros([ae, be])

    for i in range(ae):
        for j in range(be):
            if math.isnan(ed[i,j]) or sa[i,j] == 0 or scf[i,j] < 0.05: # or infr[i,j] > 0
                pass
            else:
                if ed[i,j] == 0:
                    new_ed = numpy.nanmax(ed)
                    
                    if p['sol_price'] == 0:
                        bat_price, pv_price[i,j] = costsolarbattery(new_ed, sa[i,j], scf[i,j])
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])
                    else:
                        bat_price, _ = costsolarbattery(new_ed, sa[i,j], scf[i,j])
                        pv_price[i,j] = p['sol_price']
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])
                    
                    pvInputs = [new_ed, scf[i,j], p['sol_construct_t'], 
                                pv_price[i,j], p['sol_construct_sched'], 
                                p['sol_op_years'], p['sol_op_start'], p['sol_depr_sched'],
                                p['sol_incr_cap_cost'], p['sol_fixed_om_cost'], 
                                p['sol_var_om_cost'], p['sol_decom_cost'], 
                                p['inflation_rate'], p['inflation_base_year'], 
                                p['bat_eff'], bat_price, p['bat_lifetime'], # 1, 0, 1, # 
                                p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                p['equity_rate'], p['interest_rate'], transmission_price]
                    
                    pv_price[i,j] = 0
                    lcs[i,j] = 0 # lcoe(pvInputs)
                    
                else:
                    if p['sol_price'] == 0:
                        bat_price, pv_price[i,j] = costsolarbattery(ed[i,j], sa[i,j], scf[i,j])
                    else:
                        bat_price, _ = costsolarbattery(ed[i,j], sa[i,j], scf[i,j])
                        pv_price[i,j] = p['sol_price']                    
    
                    pvInputs = [ed[i,j], scf[i,j], p['sol_construct_t'], 
                                pv_price[i,j], p['sol_construct_sched'], 
                                p['sol_op_years'], p['sol_op_start'], p['sol_depr_sched'],
                                p['sol_incr_cap_cost'], p['sol_fixed_om_cost'], 
                                p['sol_var_om_cost'], p['sol_decom_cost'], 
                                p['inflation_rate'], p['inflation_base_year'], 
                                p['bat_eff'], bat_price, p['bat_lifetime'], # 1, 0, 1, # 
                                p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                p['equity_rate'], p['interest_rate'], 0]
                    
                    lcs[i,j] = lcoe(pvInputs)

    return lcs, pv_price


# Use above functions to calculate LCOEs of wind
def parent_win(p, ed, wa, wcf, infr, trans, step):

    # Preallocate output matrices
    ae, be = numpy.shape(ed)
    lcw = numpy.zeros([ae, be])
    wind_price = numpy.zeros([ae, be])

    for i in range(ae):
        for j in range(be):
            if math.isnan(ed[i,j]) or wa[i,j] == 0 or wcf[i,j] < 0.15 or wcf[i,j] > 0.45: # or infr[i,j] > 0 # CF < 15% excluded: https://www.nrel.gov/docs/fy19osti/71814.pdf
                pass
            else: 
                if ed[i,j] == 0:
                    new_ed = numpy.nanmax(ed)

                    if p['win_price'] == 0:
                        wind_price[i,j] = costwind(new_ed, wa[i,j], wcf[i,j])
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])
                    else:
                        wind_price[i,j] = p['win_price']
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])

                    windInputs = [new_ed, wcf[i,j], p['win_construct_t'], 
                                  wind_price[i,j], p['win_construct_sched'], 
                                  p['win_op_years'], p['win_op_start'], p['win_depr_sched'],
                                  p['win_incr_cap_cost'], p['win_fixed_om_cost'], 
                                  p['win_var_om_cost'], p['win_decom_cost'], 
                                  p['inflation_rate'], p['inflation_base_year'], 
                                  1, 0, 1, # no battery inputs for wind because no storage with wind
                                  p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                  p['equity_rate'], p['interest_rate'], transmission_price]
    
                    wind_price[i,j] = 0
                    lcw[i,j] = 0 # lcoe(windInputs)
                
                else:
                    if p['win_price'] == 0:
                        wind_price[i,j] = costwind(ed[i,j], wa[i,j], wcf[i,j])
                    else:
                        wind_price[i,j] = p['win_price']
                    
                    windInputs = [ed[i,j], wcf[i,j], p['win_construct_t'], 
                                  wind_price[i,j], p['win_construct_sched'], 
                                  p['win_op_years'], p['win_op_start'], p['win_depr_sched'],
                                  p['win_incr_cap_cost'], p['win_fixed_om_cost'], 
                                  p['win_var_om_cost'], p['win_decom_cost'], 
                                  p['inflation_rate'], p['inflation_base_year'], 
                                  1, 0, 1, # no battery inputs for wind because no storage with wind
                                  p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                  p['equity_rate'], p['interest_rate'], 0]
                
                    lcw[i,j] = lcoe(windInputs)
        
    return lcw, wind_price


# Use above functions to calculate LCOEs of hydro
def parent_hyd(p, ed, flow, catcharea, slope, infr, trans, step):

    from b_calc_hydropower import calc_hydropower
   
    # Calculate total power demand (kW) from energy demand (kWh/day)
    pd = ed / 24 # kW
    new_pd = numpy.nanmean(pd)

    # Get hydro availability
    ha = calc_hydropower(catcharea, slope, flow, step) # kW/km^2
    
    # Preallocate output matrix
    ae, be = numpy.shape(ed)
    lch = numpy.zeros([ae, be])
    hydro_price = numpy.zeros([ae, be])
    hcf = numpy.zeros([ae, be])
    
    for i in range(ae):
        for j in range(be):
            if math.isnan(ed[i,j]) or ha[i,j] == 0 or slope[i,j] < 0.005: # or infr[i,j] > 0 
                # no hydraulic head < 0.6 meters (https://www.energy.gov/energysaver/planning-microhydropower-system)
                # no surface gradient < 0.5% (https://www.sciencedirect.com/science/article/pii/S1364032114003967)
                pass
            else:
                if ed[i,j] == 0:
                    new_ed = numpy.nanmax(ed)

                    # Get hydro price and capacity factor
                    if p['hyd_price'] == 0:
                        hydro_price[i,j], hcf[i,j] = costhydro(new_pd, slope[i,j], flow[i,j], ha[i,j], catcharea[i,j], step) # $/kWh
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])
                    else:
                        _, hcf[i,j] = costhydro(new_pd, slope[i,j], flow[i,j], ha[i,j], catcharea[i,j], step) # $/kWh
                        hydro_price[i,j] = p['hyd_price']
                        transmission_price = costtrans(i, j, trans, step, p['trans_price'])
                    
                    hydroInputs = [new_ed, hcf[i,j], p['hyd_construct_t'], 
                                   hydro_price[i,j], p['hyd_construct_sched'], 
                                   p['hyd_op_years'], p['hyd_op_start'], p['hyd_depr_sched'],
                                   p['hyd_incr_cap_cost'], p['hyd_fixed_om_cost'], 
                                   p['hyd_var_om_cost'], p['hyd_decom_cost'], 
                                   p['inflation_rate'], p['inflation_base_year'], 
                                   1, 0, 1, # no storage inputs for hydro because no storage with ror hydro
                                   p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                   p['equity_rate'], p['interest_rate'], transmission_price]

                    hydro_price[i,j] = 0
                    lch[i,j] = 0 # lcoe(hydroInputs)

                else:
                    # Get hydro price and capacity factor
                    if p['hyd_price'] == 0:
                        hydro_price[i,j], hcf[i,j] = costhydro(pd[i,j], slope[i,j], flow[i,j], ha[i,j], catcharea[i,j], step) # $/kWh
                    else:
                        _, hcf[i,j] = costhydro(pd[i,j], slope[i,j], flow[i,j], ha[i,j], catcharea[i,j], step) # $/kWh
                        hydro_price[i,j] = p['hyd_price']
                    
                    hydroInputs = [ed[i,j], hcf[i,j], p['hyd_construct_t'], 
                                   hydro_price[i,j], p['hyd_construct_sched'], 
                                   p['hyd_op_years'], p['hyd_op_start'], p['hyd_depr_sched'],
                                   p['hyd_incr_cap_cost'], p['hyd_fixed_om_cost'], 
                                   p['hyd_var_om_cost'], p['hyd_decom_cost'], 
                                   p['inflation_rate'], p['inflation_base_year'], 
                                   1, 0, 1, # no storage inputs for hydro because no storage with ror hydro
                                   p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                   p['equity_rate'], p['interest_rate'], 0]
                
                    lch[i,j] = lcoe(hydroInputs)
                    
                    if hydro_price[i,j] == 0:
                        lch[i,j] = 0
    
    return lch, hydro_price, hcf


# Use above functions to calculate LCOEs of wind
def parent_fossil(p, ed, infr, trans, step, ongrid_or_offgrid):

    # Preallocate output matrices
    ae, be = numpy.shape(ed)
    lcf = numpy.zeros([ae, be])

    for i in range(ae):
        for j in range(be):
            if math.isnan(ed[i,j]): # or infr[i,j] > 0
                pass
            else: 
                if ed[i,j] == 0:
                    new_ed = numpy.nanmax(ed)
                    
                    transmission_price = costtrans(i, j, trans, step, p['trans_price'])

                    # Get on- or off-grid technology inputs
                    if ongrid_or_offgrid == "on":
                        Inputs = [new_ed, p['col_cf'], p['col_construct_t'], 
                                      p['col_price'], p['col_construct_sched'], 
                                      p['col_op_years'], p['col_op_start'], p['col_depr_sched'],
                                      p['col_incr_cap_cost'], p['col_fixed_om_cost'], 
                                      p['col_var_om_cost'], p['col_decom_cost'], 
                                      p['inflation_rate'], p['inflation_base_year'], 
                                      1, 0, 1, # no battery inputs
                                      p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                      p['equity_rate'], p['interest_rate'], transmission_price]
                    else:
                        Inputs = [new_ed, p['dies_cf'], p['dies_construct_t'], 
                                      p['dies_price'], p['dies_construct_sched'], 
                                      p['dies_op_years'], p['dies_op_start'], p['dies_depr_sched'],
                                      p['dies_incr_cap_cost'], p['dies_fixed_om_cost'], 
                                      p['dies_var_om_cost'], p['diesl_decom_cost'], 
                                      p['inflation_rate'], p['inflation_base_year'], 
                                      1, 0, 1, # no battery inputs
                                      p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                      p['equity_rate'], p['interest_rate'], transmission_price]

                    lcf[i,j] = 0 # lcoe(Inputs)
                
                else:
                    # Get on- or off-grid technology inputs
                    if ongrid_or_offgrid == "on":
                        Inputs = [ed[i,j], p['col_cf'], p['col_construct_t'], 
                                      p['col_price'], p['col_construct_sched'], 
                                      p['col_op_years'], p['col_op_start'], p['col_depr_sched'],
                                      p['col_incr_cap_cost'], p['col_fixed_om_cost'], 
                                      p['col_var_om_cost'], p['col_decom_cost'], 
                                      p['inflation_rate'], p['inflation_base_year'], 
                                      1, 0, 1, # no battery inputs
                                      p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                      p['equity_rate'], p['interest_rate'], 0]
                    elif ongrid_or_offgrid == "off":
                        Inputs = [ed[i,j], p['dies_cf'], p['dies_construct_t'], 
                                      p['dies_price'], p['dies_construct_sched'], 
                                      p['dies_op_years'], p['dies_op_start'], p['dies_depr_sched'],
                                      p['dies_incr_cap_cost'], p['dies_fixed_om_cost'], 
                                      p['dies_var_om_cost'], p['dies_decom_cost'], 
                                      p['inflation_rate'], p['inflation_base_year'], 
                                      1, 0, 1, # no battery inputs
                                      p['tax_rate'], p['debt_frac'], p['debt_rate'], 
                                      p['equity_rate'], p['interest_rate'], 0]

                    lcf[i,j] = lcoe(Inputs)
        
    return lcf



