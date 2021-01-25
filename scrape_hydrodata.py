# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:32:45 2020

Scrape hydropower plant locations data from Guatemalan govt website.

@author: clhenry
"""

from bs4 import BeautifulSoup
import requests
import pandas as pd

# Call url and collect html text
response = requests.get("http://www.cnee.gob.gt/wp/?page_id=239")
soup = BeautifulSoup(response.text)

# Get all href links on page
links = []
for link in soup.find_all('a'):
    links.append(link.get('href'))

# Clean list of links
links = [x for x in links if x != None]

# Get indices of relevant links (links to dam pages)
index_start = [idx for idx, s in enumerate(links) if 'maps' in s][0] + 1
index_end = [idx for idx, s in enumerate(links) if '#' in s][0]

# Create new list of links to only relevant pages
pages = links[index_start:index_end]

# Create empty dictionary
df = []

# Iterate through each sub-link
for page in pages:
    # Call urls and collect html text
    page_response = requests.get(page)
    page_soup = BeautifulSoup(page_response.text)

    # Get all table text on page
    txts = []
    for i in range(len(page_soup.find_all('td'))):
        txts.append(page_soup.find_all('td')[i].get_text())

    # Collect relevant info into table
    name_index = [idx for idx, s in enumerate(txts) if 'presa' in s.lower()][1] + 1
    lat_index = [idx for idx, s in enumerate(txts) if 'latitud' in s.lower()][0] + 1
    lon_index = [idx for idx, s in enumerate(txts) if 'longitud' in s.lower()][0] + 1
    flow_index = [idx for idx, s in enumerate(txts) if 'caudal a turbinar' in s.lower()][0] + 1
    turbine_index = [idx for idx, s in enumerate(txts) if 'tipo de turbina' in s.lower()][0] + 1
    nameplate_index = [idx for idx, s in enumerate(txts) if 'potencia efectiva' in s.lower()][0] + 1
    
    df.append({
            'Name': txts[name_index],
            'Latitude': txts[lat_index],
            'Longitude': txts[lon_index],
            'Turbine Flow': txts[flow_index],
            'Turbine Type': txts[turbine_index],
            'Capacity': txts[nameplate_index]
        })

# Convert dictionaries into dataframe
dictionary = pd.DataFrame(df)

# Write as csv
dictionary.to_excel('guatemala_hydroplants.xlsx', index = False, header=True)
