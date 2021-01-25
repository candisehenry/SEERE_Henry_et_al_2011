# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:12:06 2020

Set up data inputs for SEERE run.
Hardcoded for Guatemala.

@author: clhenry
"""

def d_setup(run_d_setup, flow_suffix, temp_suffix, wind_multiplier):

    from set_coordinates import set_coordinates
    from a_get_hydrodata import get_hydrodata
    from b_calc_solarpower import calc_solar
    from b_calc_windpower import calc_wind_tif_v2
    from b_calc_enedem import calc_energydemand_tif
    from c_interp_solar import interp_solar
    from c_infrastructure import get_infrastructure
#   from b_calc_windpower import calc_wind
#   from b_calc_windpower import calc_wind_tif
#   from b_calc_enedem import calc_energydemand_csv
#   from c_interp_wind import interp_wind_csv
#   from c_interp_enedem import interp_enedem

    # Size of grid cells
    step = 1/120 # step size we want
    step_wind = 0.00417
    step_hyd = 0.00834
    step_pop = 0.00834 # gtm input
    step_infr = 0.004167
    step_trans = 0.004167

    # Set coordinates (get coordinates from qgis)
    ## standardize all sets to solar coordinates
    coord = set_coordinates('Guatemala', 'POWER_Regional_Interannual_199001_201912_c2f60e61.csv')
    coord_win = [13.7382822040000008, 17.818714142000001, -92.2223587040000012, -88.2256698610000001]
    coord_hyd = [13.7382822040000008, 17.8248822040000015, -92.2223587040000012, -88.2191587039999945] # [13.7337510510000005,15.2683343229999995,-92.2345844440000064,-89.2012512320000042]
    coord_pop = [13.7375166666666679, 17.8208500000000001, -92.2250700000000023, -88.2250700000000023] # [13.7375109900000023, 17.8208475900000032, -92.2250737800000309, -88.2250705800000219] # [13.7412501466431500, 17.8245834636431510, -92.2262496413472235, -88.2345829906472261]  
    coord_infr = [13.6924346630189255, 17.8674380030189255, -92.3816732455897096, -88.0525031155897153]
    coord_trans = [13.738644377, 17.817811043999999, -92.2293775180000068, -88.2335441849999995]


    # Run all calculations and interpolations
    if run_d_setup == 'n' or run_d_setup == 'N':
        pass

    else:
        # Run solar calculations and interpolations
        calc_solar(folder = 'Guatemala', infile = 'POWER_Regional_Interannual_199001_201912_c2f60e61.csv', 
                   outfile1 = 'guatemala_pv_power.txt', outfile2 = 'guatemala_solar_cf.txt', run = temp_suffix)
        interp_solar(folder = 'Guatemala', infile = 'guatemala_pv_power.txt', 
                     outfile = 'guatemala_pv_interp.npy', step = step, coord = coord, run = temp_suffix)
        interp_solar(folder = 'Guatemala', infile = 'guatemala_solar_cf.txt', 
                     outfile = 'guatemala_scf_interp.npy', step = step, coord = coord, run = temp_suffix)
        
        # Run wind calculations (interpolated)
        calc_wind_tif_v2(folder = 'Guatemala', infile_wndspd = 'GlobalAtlasWind_GTM_windspeed_50m.tif', 
                        infile_temp = 'gridded_mean_temp.npy',
                        outfile1 = 'guatemala_wind_interp.npy', outfile2 = 'guatemala_wcf_interp.npy', 
                        step_from_wndspd = step_wind, step_from_temp = step,
                        step_to = step, 
                        coord_from_wndspd = coord_win, coord_from_temp = coord,
                        coord_to = coord,
                        variability = wind_multiplier)
        
        # Run energy demand calculations for solar and wind
        calc_energydemand_tif(folder = 'Guatemala', infile = 'hogares_offgrid_munis.tif', 
                            outfile = 'guatemala_demand_interp_survey.npy', 
                            step_from = step_pop, step_to = step, 
                            coord_from = coord_pop, coord_to = coord)
        
        
        
        # Run hydro calculations and interpolations (should already be in correct step size)
        get_hydrodata(folder = 'Guatemala', catcharea_infile = 'guatemala_catcharea.tif', 
                      slope_infile = 'guatemala_slope.tif', flow_infile = 'guatemala_' + flow_suffix + '_flow.tif', 
                      catcharea_outfile = 'catcharea.npy', slope_outfile = 'slope.npy', 
                      flow_outfile = 'flow.npy')
        
        # # Run energy demand calculations for hydro
        # calc_energydemand_tif(folder = 'Guatemala', infile = 'hogares_offgrid_munis.tif', 
        #                     outfile = 'guatemala_demand_interp_survey.npy', 
        #                     step_from = step_pop, step_to = step_hyd, 
        #                     coord_from = coord_pop, coord_to = coord_hyd)



        # Run infrastructure set up
        get_infrastructure(folder = ['Guatemala', 'Guatemala/LandCover/'], 
                           infile = ['urbanareas_new.tif'], 
                           outfile = 'guatemala_urbanareas.npy',
                           step_from = [step_infr], step_to = step_hyd,
                           coord_from = [coord_infr], coord_to = coord_hyd)

        get_infrastructure(folder = ['Guatemala', 'Guatemala/Transmission/'], 
                           infile = ['transmission.tif'], 
                           outfile = 'guatemala_transmission.npy',
                           step_from = [step_trans], step_to = step_hyd,
                           coord_from = [coord_trans], coord_to = coord_hyd)

#       calc_wind(folder = 'Guatemala', infile = 'POWER_Regional_Interannual_199001_201912_c2f60e61.csv', 
#                 outfile1 = 'guatemala_wind_power.txt', outfile2 = 'guatemala_wind_cf.txt')
#       interp_wind_csv(folder = 'Guatemala', infile = 'guatemala_wind_power.txt', 
#                   outfile = 'guatemala_wind_interp.npy', step = step, coord = coord)
#       interp_wind_csv(folder = 'Guatemala', infile = 'guatemala_wind_cf.txt', 
#                  outfile = 'guatemala_wcf_interp.npy', step = step, coord = coord)
#
#       # Run wind calculations (interpolated )
#       calc_wind_tif(folder = 'Guatemala', infile_wndspd = 'GlobalAtlasWind_GTM_windspeed_50m.tif', 
#                  infile_pwrdns = 'GlobalAtlasWind_GTM_powerdens_50m.tif',
#                  outfile1 = 'guatemala_wind_interp.npy', outfile2 = 'guatemala_wcf_interp.npy', 
#                  step_from = step_wind, step_to = step, 
#                  coord_from = coord_win, coord_to = coord,
#                  variability = wind_multiplier)
#
#       calc_energydemand_csv(folder = 'Guatemala', infile = 'popdat.csv', 
#                          outfile = 'guatemala_demand.txt')
#       interp_enedem(folder = 'Guatemala', infile = 'guatemala_demand.txt', 
#                          outfile = 'guatemala_demand_interp.npy', 
#                          step = step, coord = coord)
#
#       # Run energy demand calculations and interpolations using high-resolution tiff file
#       calc_energydemand_npy(folder = 'Guatemala', infile = 'gtm_ppp_2020_1km_EDIT.tif', 
#                          outfile = 'guatemala_demand_interp_gtm.npy', 
#                          step_from = step_pop, step_to = step, 
#                          coord_from = coord_pop, coord_to = coord)
#
#       calc_energydemand_npy(folder = 'Guatemala', infile = 'gtm_ppp_2020_1km_EDIT.tif', 
#                          outfile = 'guatemala_demand_hydro_interp_gtm.npy', 
#                          step_from = step_pop, step_to = step, 
#                          coord_from = coord_pop, coord_to = coord_hyd)


    return step_hyd, coord, coord_hyd
    # return step, coord, coord_hyd
