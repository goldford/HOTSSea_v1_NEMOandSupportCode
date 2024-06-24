''' convert DFO CUR data
    G Oldford Jul 2022
    
    Based on convert_HakaiCTDs.py
    data in: large NC file of DFO IOS Moored current metre
    data out: many .NC files (one per cast)

    Source data from 1978 to 2022 downloaded from CIOOS Pacific
    https://data.cioospacific.ca/erddap/tabledap/IOS_CUR_Moorings.html
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

indir = '/project/6006412/goldford/data/evaluation/IOS_MooredCUR_Data/'
savedir = '/project/6006412/goldford/data/evaluation/IOS_MooredCUR_Data/prepped/'
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
        
        time = nc.num2date(data1['time'][:], data1['time'].units)
        
        T = data1['TEMPPR01'][:].flatten();
        T[T == 99999.] = np.nan
        T[T == -99999.] = np.nan
        T[T == -99.] = np.nan
        
        P = data1['PRESPR01'][:].flatten();
        P[P == 99999.] = np.nan
        P[P == -99999.] = np.nan
        P[P == -99.] = np.nan
        
        depth = data1['depth'][:].flatten();
        depth[depth == 99999.] = np.nan
        depth[depth == -99999.] = np.nan
        depth[depth == -99.] = np.nan
        
        # eastward vel
        u = data1['LCEWEL01'][:].flatten();
        u[u == 99999.] = np.nan
        u[u == -99999.] = np.nan
        u[u == -99.] = np.nan
        u[u < -99900.] = np.nan
        u[u > 99900.] = np.nan
        
        # northward vel
        v = data1['LCNSEL01'][:].flatten();
        v[v == 99999.] = np.nan
        v[v == -99999.] = np.nan
        v[v == -99.] = np.nan
        v[v < -99900.] = np.nan
        v[v > 99900.] = np.nan
        
        # upward vel
        z = data1['LRZASP01'][:].flatten();
        z[z == 99999.] = np.nan
        z[z == -99999.] = np.nan
        z[z == -99.] = np.nan
        z[z < -99900.] = np.nan
        z[z > 99900.] = np.nan
        
          
        for m_id in unique_mooringids:
          
          mask = (data1['profile'][:] == m_id)

          T_m = T[mask]
          P_m = P[mask]
          u_m = u[mask]
          v_m = v[mask]
          z_m = z[mask]
          
          # only need one val for these 
          lon_m = lon[mask]
          lat_m = lat[mask]
          dep_m = depth[mask]
          lat_m1 = lat_m[0]
          lon_m1 = lon_m[0]
          dep_m1 = dep_m[0]
          
          projects_m = projects[mask]
          instDep_m = instDep[mask]
          instDep_m1 = instDep_m[0]

          time_m = time[mask]
          timeRange = [time_m[0], time_m[-1]]
          # date range in seconds
          dts = np.array([(DT - time_m[0]).total_seconds() for DT in time_m])
          dt = np.mean(np.diff(dts))

          prj = str(projects_m[0])
          prj_id = prj + str(m_id)
          
          # replace special chars
          replaced_prj = re.sub('[,. &\-/\\\\]', '_', prj_id)
          
          print("writing to: ", replaced_prj)
          
          #CurrentMeter(self, saveDir, station, timeRange, instDep, dts, T, P, u, v, w, dt, dep, lat, lon):
          
          outfile = write_to_netcdf().CurrentMeter(savedir, replaced_prj, timeRange, instDep_m1, dts, T_m, P_m, u_m, v_m, z_m, dt, dep_m1, lat_m1, lon_m1)
          print("Saved", outfile)
