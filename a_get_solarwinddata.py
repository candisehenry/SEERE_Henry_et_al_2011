# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:22:47 2020

Functions to import NASA solar and wind data.

@author: clhenry
"""

import requests
import pandas as pd
import csv


# Download NASA data using API
def import_nasa_api(dataURL):

    """
    Documentation: https://power.larc.nasa.gov/docs/services/api/v1/temporal/climatology/
    Get coordinates: https://power.larc.nasa.gov/data-access-viewer/

    More info on APIs in Python: https://rapidapi.com/blog/how-to-use-an-api-with-python/
    """
    
    # URL for Guatemala monthly min/max/mean temp/wind speed/insolation data for 10 years
    if dataURL == '':
        dataURL = 'https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?&request=execute&identifier=Regional&parameters=T2M_MAX,T2M_MIN,T2M,ALLSKY_SFC_SW_DWN,WS50M_MIN,WS50M_MAX,WS50M&startDate=1990&endDate=2019&userCommunity=SSE&tempAverage=INTERANNUAL&bbox=13.50,-92.50,18.00,-88.00&outputList=CSV'
    
    # Make sure grid coordinates are larger than country extent because values will be rounded to nearest 25 degrees
    
    # Get API
    response = requests.get(dataURL)
    
    # Check response
    print(response) # 200 means successful
    
    # Output csv file with data
    results = response.json()
    
    return results['outputs']
    # Then go to output csv link in browser (to download file)


# Read in csv file of raw NASA data
def import_raw_nasa_data(file):
    f = open(file)
    rdr = csv.reader(f)
    
    data = []
    
    #Throw away all lines up to and include the line that has '-END HEADER-' in the first cell of the line
    while True:
        line = next(rdr)
        if line[0] == '-END HEADER-':
            break
        
    # Now take the header row
    line = next(rdr)
    data.append(line)

    # Now take all non-blank lines
    while True:
        try:
            line = next(rdr)
            if line[0] != '':
                data.append(line)
        except:
            break

    data_array = pd.DataFrame(data[1:], columns = data[0], dtype = float)

    return data_array