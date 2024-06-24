''' convert DFO CTD data
    G Oldford Jul 2022
    
    Based on convert_HakaiCTDs.py
    data in: large NC file of DFO IOS Moored CTD casts
    data out: many .NC files (one per cast)

    Source data from 1978 to 2022 downloaded from CIOOS Pacific
    https://data.cioospacific.ca/erddap/tabledap/IOS_CTD_Moorings.html
    Used Salish Sea bounding box but this included stations off WCVI
    
'''

import glob
import os
import gsw
import numpy as np
import pandas as pd
import analysispkg
import re
import netCDF4 as nc

import sys
sys.path.insert(1, '/project/6006412/goldford/py3-8-10_Jul9/lib/python3.8/site-packages/python-analysis-package/etc/data_conversion/')
from data_conversion_functions import write_to_netcdf

indir = '/project/6006412/goldford/data/evaluation/IOS_Moored_Data/'
savedir = '/project/6006412/goldford/data/evaluation/IOS_Moored_Data/MCTD_prepped/'
os.makedirs(savedir, exist_ok=True)
files = sorted(glob.glob(os.path.join(indir, "*.nc")))

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



# main loop 
for file in files: 
    print("Working on ", file)
    with nc.Dataset(file, 'r') as data1:
        
        mooring_id = data1['profile'][:].flatten();
        projects = data1['project'][:].flatten();
        unique_mooringids = np.unique(mooring_id)
        
        lat = data1['latitude'][:].flatten();
        lon = data1['longitude'][:].flatten();
        instDep = data1['instrument_depth'][:].flatten();
        dep = -999.
        
        time = nc.num2date(data1['time'][:], data1['time'].units)
        
         # the columns for temp, salin, etc are not consistently populated
        T = data1['sea_water_temperature'][:].flatten();
        T[T == 99999.] = np.nan
        T2 = data1['TEMPST01'][:].flatten();
        T2[T2 == 99999.] = np.nan
        T3 = data1['TEMPS901'][:].flatten();
        T3[T3 == 99999.] = np.nan
        T4 = data1['TEMPS902'][:].flatten();
        T4[T4 == 99999.] = np.nan
        T5 = data1['TEMPS601'][:].flatten();
        T5[T5 == 99999.] = np.nan
        T6 = data1['TEMPS602'][:].flatten();
        T6[T6 == 99999.] = np.nan
        
        # ////////////////////
        # replace nans with alt vals if avail
        T[np.isnan(T)] = T2[np.isnan(T)]
        T[np.isnan(T)] = T3[np.isnan(T)]
        T[np.isnan(T)] = T4[np.isnan(T)]
        T[np.isnan(T)] = T5[np.isnan(T)]
        T[np.isnan(T)] = T6[np.isnan(T)]
        
        S = data1['sea_water_practical_salinity'][:].flatten();
        S[S == 99999.] = np.nan
        S1 = data1['PSALST01'][:].flatten();
        S1[S1 == 99999.] = np.nan       
        S2 = data1['PSALST02'][:].flatten();
        S2[S2 == 99999.] = np.nan
        S[np.isnan(S)] = S1[np.isnan(S)]
        S[np.isnan(S)] = S2[np.isnan(S)]
        
        P = data1['sea_water_pressure'][:].flatten();
        P[P == 99999.] = np.nan
        P2 = data1['PRESPR01'][:].flatten();
        P2[P2 == 99999.] = np.nan  
        P[np.isnan(P)] = P2[np.isnan(P)]
          
        for m_id in unique_mooringids:
          
          mask = (data1['profile'][:] == m_id)

          T_m = T[mask]
          S_m = S[mask]
          P_m = P[mask]
          lon_m = lon[mask]
          lat_m = lat[mask]
          lat_m1 = lat_m[0]
          lon_m1 = lon_m[0]
          
          projects_m = projects[mask]
          instDep_m = instDep[mask]
          instDep_m1 = instDep_m[0]

          absSal = gsw.conversions.SA_from_SP(S_m, P_m, lon_m, lat_m)
          consTemp = gsw.conversions.CT_from_t(absSal, T_m, P_m)
          pDen = gsw.rho(absSal, consTemp, np.zeros_like(P_m))

          time_m = time[mask]
          timeRange = [time_m[0], time_m[-1]]
          dts = np.array([(DT - time_m[0]).total_seconds() for DT in time_m])
          dt = np.mean(np.diff(dts))

          prj = str(projects_m[0])
          prj_id = prj + str(m_id)
          
          # replace special chars
          replaced_prj = re.sub('[,. &-]', '_', prj_id)
          
          print("writing to: ", replaced_prj)
          outfile = write_to_netcdf().MooredCTD(savedir, replaced_prj, timeRange, instDep_m1, dts, T_m, S_m, P_m, consTemp, absSal, pDen, dt, dep, lat_m1, lon_m1)
          print("Saved", outfile)
