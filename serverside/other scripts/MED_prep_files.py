''' convert DFO IOS Waterproperties MED files
    G Oldford Jul 2022
    
    Based on code in python-analysis-package etc/data_conversion_functions.py IOS.CurrentMeter
    data in: many MED files from the waterproperties.ca website (see metadata)
    data out: many .NC files (one per MED measurement)
    
    JUST STARTED - PUTTING ON HOLD
    
'''

import glob
import os
import gsw
import numpy as np
import pandas as pd
import analysispkg
import re
import netCDF4 as nc
import datetime as datetime

import sys
sys.path.insert(1, '/project/6006412/goldford/py3-8-10_Jul9/lib/python3.8/site-packages/python-analysis-package/etc/data_conversion/')
from data_conversion_functions import write_to_netcdf

indir = '/project/6006412/goldford/data/evaluation/WaterProp_MED/'
savedir = '/project/6006412/goldford/data/evaluation/WaterProp_MED/prepped/'
os.makedirs(savedir, exist_ok=True)
files = sorted(glob.glob(os.path.join(indir, "*.MED")))

n=0
badvals=0

# helper to filter profiles with possible null values
# where all values are the same (e.g. -99.847)
def is_unique(s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()

# replaces bad values like -9.99e-29
def replace_bad(bad_val, T, SP, P, badvals):
    ikeep = ~ (np.isclose(SP, bad_val, bad_val*2) | np.isclose(T, bad_val, bad_val*2))
    badvals += np.sum(~ikeep)
    T,SP,P=T[ikeep],SP[ikeep],P[ikeep]
    return T,SP,P,badvals

i = 0

# main loop 
for file in files: 
    print("Working on ", file)
    
    with open(file,'r',encoding='latin-1') as f:
      lines = [line for line in f]
    
    timeRange = []    
    
    # lN = line number
    for lN,line in enumerate(lines):
      
      
      if len(timeRange) < 2:
        if 'START TIME' in line:
          #sT = line[line.find('UTC') + 4:line.rfind('.')]
          sT = line[line.find('GMT') + 4:line.rfind('.')]
          sT = datetime.datetime.strptime(sT,'%Y/%m/%d %H:%M:%S')
          # GMT 1979/08/19 00:00:00.000
          timeRange.append(sT)
          
          print(sT)
          
        if 'END TIME' in line:
          #eT = line[line.find('UTC') + 4:line.rfind('.')]
          eT = line[line.find('GMT') + 4:line.rfind('.')]
          eT = datetime.datetime.strptime(eT,'%Y/%m/%d %H:%M:%S')
          timeRange.append(eT)
          
          print(eT)
      
    # for debug  
    #print(lN, line)
    i = i + 1
    if i == 1:
      break
    

          
    #outfile = write_to_netcdf().CurrentMeter(savedir, replaced_prj, timeRange, instDep_m1, dts, T_m, P_m, u_m, v_m, z_m, dt, dep_m1, lat_m1, lon_m1)
    #print("Saved", outfile)
