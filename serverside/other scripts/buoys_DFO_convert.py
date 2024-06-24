''' convert DFO buoy data processed from original ECCC
    G Oldford Jun 8 2022

    data in: .csv file of all buoy data in in ISO 8859-1 (two header rows, second row = units)
    data out: many .NC files (one per buoy)

    NetCDF data downloaded by G Oldford
    from https://catalogue.cioospacific.ca/dataset/ca-cioos_b9c71eb2-b750-43d5-a50a-aee173916736
    
    same data as convert_BUOYS_ECCC.py processes except this is in CSV w/ different data attributes
    
    To do: 
      probably better not to read data with pandas (slow)
      instead try approach used in convert_BUOYS_ECCC.py
    
'''

import glob
import os
import gsw
import numpy as np
import pandas as pd
import analysispkg
import re
import math
import datetime

import sys
sys.path.insert(1, '/project/6006412/goldford/py3-8-10_Jul9/lib/python3.8/site-packages/python-analysis-package/etc/data_conversion/')
from data_conversion_functions import write_to_netcdf

indir = '/project/6006412/goldford/data/evaluation/Buoys/'
savedir = '/project/6006412/goldford/data/evaluation/Buoys/prepped/'

os.makedirs(savedir, exist_ok=True)
files = sorted(glob.glob(os.path.join(indir, "*.csv")))

      
# main loop 
for file in files: 
    print("Working on ", file)
    d = pd.read_csv(file, skiprows=range(1,2), dtype={"STN_ID": object, 
                                                      "latitude": float, "longitude": float,
                                                      "water_depth": float, "VCAR": float,
                                                      "VTPK": float, "VWH": float, 
                                                      "VCMX": float, "VTP": float,
                                                      "WDIR": float, "WSPD": float, 
                                                      "WSS": float, "GSPD": float,
                                                      "WDIR2": float, "WSPD2": float, "WSS_2": float, 
                                                      "GSPD2": float, "ATMS": float, "ATMS2": float, 
                                                      "DRYT": float, "SSTP": float, "SSTP_flags": float
                                                      },
                                                      low_memory=False)
    
    print("Reading Lighthouse data, length: " + str(len(d)))
    print(d.columns)
    buoys = d[['STN_ID']].drop_duplicates()
    
    #d['date'] = pd.to_datetime(d['time'])
    #d['year'] = d['date'].dt.year
    
    print("Unique buoys: ", buoys)
    print("Length of file: ", len(d))
    
    # first file row contains vars, units in second row
    d2 = pd.read_csv(file, nrows=2)
    d_units = d2.iloc[0]
    print("Units: ", d_units)
    
    
    # for each buoy
    for idx, row in buoys.iterrows():
        
        buoy_ID = row['STN_ID']
        buoy_data = d[(d['STN_ID'].values == buoy_ID)]
        print("Working on: ", buoy_ID)
        print("length of single buoy data rec: ", len(buoy_data))
        
        lats = []; lons = []; Ts = []; dates1 = []
        
        # for each row in buoy dataset
        for idx2, row2 in buoy_data.iterrows():
          
          
          lat = row2['latitude']
          lon = row2['longitude']
          #print("lat ", lat)
          #print("lon ", lon)
        
          T = row2['SSTP']
          # error catching, skipping records if T is bad
          try:
           if T is None: 
             #print('T is None')
             continue
             
          except NameError:
             #print ("T is not defined")
             continue

          if math.isnan(T):
            #print ("T is nan")
            continue
            
          if T < 0:
            #print ("T is < 0")
            continue
            
          Ts.append(T); #lats.append(lat); lons.append(lon) - assuming lats and lons don't change over time
          dates1.append (datetime.datetime.strptime(row2['time'],'%Y-%m-%dT%H:%M:%SZ') )

        # Write the file
        if len(Ts) > 0:
        
          # Station name - no 'name' in the DFO version of ECCC buoy data, just code
          #stn_name = stn_name.replace(' ' , '-')
          span = dates1[0].strftime('%Y-%m') + '_' + dates1[-1].strftime('%Y-%m')
          #outfile = os.path.join(dir_dest, '_'.join([stn_id, stn_name, span]) + '.nc' )
          outfile = os.path.join(savedir, '_'.join([buoy_ID, span]) + '.nc')
        
          print("lengthchecks: ")
          print(len(dates1))
          print(len(Ts))
        
          print("outfile: ", outfile)
          write_to_netcdf().SST(outfile, lon, lat, dates1, Ts, 'degC')
        
        else: 
          print("no sst records found, continuing")
          continue
        



        
        



      