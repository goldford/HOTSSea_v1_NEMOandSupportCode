''' convert DFO CTD data
    G Oldford Jul 2021

    Based on convert_HakaiCTDs.py and other FA12 code
    files in: CSV of DFO IOS CTD rosette bottle samples
    files out: .NC files (one per cast, separate files for temp salin)

    Data from 1978 to 2020 downloaded from CIOOS Pacific in July 2021
    https://data.cioospacific.ca/erddap/tabledap/IOS_BOT_Profiles.html
    Used Salish Sea bounding box but this included stations off WCVI

    Notes: Rosette data includes nitrate / nitrite, oxygen conc, chlorophyll-a, 
    sea water conductivity, phosphate, salinity, temperature, silicate
    but these data are not currently written out
    
    Salinity and temp profiles are inconsistent and depends on equip and era. 
    
'''

# to do: Aug 17 2021 write NC files with conservative 
#        temperature when salinity and temp are both present
import netCDF4 as nc
import glob
import os
from os.path import abspath, expanduser, exists, join
import gsw
import numpy as np
import pandas as pd
import analysispkg
import datetime as datetime

import sys
sys.path.insert(1, '/project/6006412/goldford/py3-8-10_Jul9/lib/python3.8/site-packages/python-analysis-package/etc/data_conversion/')
from data_conversion_functions import write_to_netcdf

indir = '/project/6006412/goldford/data/evaluation/IOS_Rosette_CTD/'
savedir = '/project/6006412/goldford/data/evaluation/IOS_Rosette_CTD/prepped/'

os.makedirs(savedir, exist_ok=True)

files = sorted(glob.glob(os.path.join(indir, "*.csv")))
n=0
n_temp=0
n_salin=0
stationnodata_T=0
stationnodata_S=0
stationnodata_P=0
nullvals_T=0
nullvals_S=0
nbadvals_t=0
nbadvals_s=0

# helper to filter profiles with possible null values
# where all values are the same (e.g. -99.847)
def is_unique(s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()

# drop null values (helper by GO)
def drop_nulls(V, P, datesnp, badvals):
    ikeep = ~ (np.isnan(V))
    badvals += np.sum(~ikeep)
    V=V[ikeep]
    P=P[ikeep]
    datesnp=datesnp[ikeep]
    return V, P, datesnp, badvals
    
# replaces nodata values like -9.99e-29
def drop_bad(bad_val, V, P, datesnp, nbadvals):
    ikeep = ~ (np.isclose(V, bad_val, bad_val*2))
    nbadvals += np.sum(~ikeep)
    V=V[ikeep]
    P=P[ikeep]
    datesnp=datesnp[ikeep]
    return V, P, datesnp, nbadvals
    
# helper from data_conversion_functions.py
# used by function above - should revert if code is integrated 
def makePath(path):
    return abspath(expanduser(path)) 
  
# based on convert_Hakai_CTDs
def Rosette_Bottle_Temp(saveDir, station, time, T, D, lat, lon):
    writeName = station.replace('/', '-') + '_Temp' + '.nc'
    
    # function duplicated here for convenience - should be switched back
    #svPath = utils().makePath(join(saveDir, writeName))
    svPath = makePath(join(saveDir, writeName))
    with nc.Dataset(svPath, 'w', format='NETCDF4') as f:
        # dimensions
        f.createDimension('single', 1)
        f.createDimension('time_counter', None)

        # variables
        # time_counter
        time_counter = f.createVariable('time', 'float32', ('time_counter'))
        time_counter.units = 'seconds since 1900-01-01 00:00:00'
        time_counter[:] = nc.date2num(time, time_counter.units)

        Temp = f.createVariable('temperature', 'f8', 'time_counter')
        Temp[:] = T
        Temp.type = 'in-situ'
        Temp.units = 'degC'
        Temp.description = 'recorded temperature'

       # cT = f.createVariable('cTemp', 'f8', 'time_counter')
       # cT[:] = consTemp
       # cT.units = 'degC'
       # cT.description = 'conservative temperature'
       
        depth = f.createVariable('depth', 'f8', 'time_counter')
        depth[:] = D
        depth.units = 'm'
        depth.description = 'water depth at cast location'

        Lat = f.createVariable('latitude', 'f4', 'single')
        Lat[:] = lat
        Lat.units = 'degrees north'
        Lat.description = 'sample latitude'

        Lon = f.createVariable('longitude', 'f4', 'single')
        Lon[:] = lon
        Lon.units = 'degrees east'
        Lon.description = 'sample longitude (-180 -> 180 degrees)'

        return svPath

def Rosette_Bottle_Salin(saveDir, station, time, S, D, lat, lon):
    writeName = station.replace('/', '-') + '_Salin' + '.nc'

    # function code duplicated above for debug - should be switched back
    #svPath = utils().makePath(join(saveDir, writeName))
    svPath = makePath(join(saveDir, writeName))

    with nc.Dataset(svPath, 'w', format='NETCDF4') as f:
        # dimensions
        f.createDimension('single', 1)
        f.createDimension('time_counter', None)

        # variables
        # time_counter
        time_counter = f.createVariable('time', 'float32', ('time_counter'))
        time_counter.units = 'seconds since 1900-01-01 00:00:00'
        time_counter[:] = nc.date2num(time, time_counter.units)

        Sal = f.createVariable('salinity', 'f8', 'time_counter')
        Sal[:] = S
        Sal.type = 'practical'
        Sal.units = 'psu'
        Sal.description = 'recorded salinity'

       # aSal = f.createVariable('aSal', 'f8', 'time_counter')
       # aSal[:] = absSal
       # aSal.units = 'g/kg'
       # aSal.description = 'absolute salinity'

        d = f.createVariable('depth', 'f8', 'time_counter')
        d[:] = D
        d.units = 'm'
        d.description = 'water depth at cast location'

        Lat = f.createVariable('latitude', 'f4', 'single')
        Lat[:] = lat
        Lat.units = 'degrees north'
        Lat.description = 'sample latitude'

        Lon = f.createVariable('longitude', 'f4', 'single')
        Lon[:] = lon
        Lon.units = 'degrees east'
        Lon.description = 'sample longitude (-180 -> 180 degrees)'

        return svPath

# main loop (just one file r now)
for file in files: 
    print("Working on ", file)
    d = pd.read_csv(file, skiprows=range(1,2), dtype={"profile": object, "instrument_model": object, 
                                                      "instrument_serial_number": object, 
                                                      "latitude": float, "longitude": float,
                                                      "geographic_area": object,
                                                      "depth": float, "CNDCST01": float,
                                                      "DOXMZZ01": float, "DOXYZZ01": float,
                                                      "NTRZAAZ1": float, "PHOSAAZ1": float,
                                                      "PRESPR01": float, "PSALST01": float,
                                                      "PSALST02": float, "PSALBST01": float,
                                                      "sea_water_practical_salinity": float, 
                                                      "sea_water_pressure": float, "sea_water_temperature": float,
                                                      "SSALBST01": float, "SSALST01": float, 
                                                      "SLCAAAZ1": float,
                                                      "TEMPRTN1": float, 
                                                      "TEMPS601": float, "TEMPS602": float, 
                                                      "TEMPS901": float, "TEMPS902": float, 
                                                      "TEMPST01": float
                                                      })
    print(d.columns)
    profiles = d[['profile']].drop_duplicates()
    print("Unique casts: ", len(profiles))
    print(d.head())
   
    # units are in second row
    d_units = pd.read_csv(file, nrows=1)
    print("Units: ", d_units)
    
    # temp
    for idx, row in profiles.iterrows():
        
        castpk = row['profile']
        p = d[d['profile'].values == castpk]
        dep = 0
        
        bottle_str = str(p['profile'].values[0])
        lat = p['latitude'].values[0]
        lon = p['longitude'].values[0]
         
        # depth / pressure
        # sea_water_pressure (dbar), PRESPR01 (decibar), depth (m)
        if (p['sea_water_pressure'].isnull().sum() < (0.5*len(p))) & (is_unique(p['sea_water_pressure']) == False):
            P = p['sea_water_pressure'].values
        elif (p['PRESPR01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['PRESPR01']) == False):
            P = p['PRESPR01'].values
        else:
            stationnodata_P += 1
            continue
         
        # temperature - check that more than half are populated and not the same value, convert
        # sea_water_temperature (degC), TEMPS601 (degC), TEMPS602 (degC), TEMPS901 (degC), TEMPS902 (degC), TEMPST01 (degC)
        if (p['sea_water_temperature'].isnull().sum() < (0.5*len(p))) & (is_unique(p['sea_water_temperature']) == False):
            T = p['sea_water_temperature'].values
        elif (p['TEMPS601'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPS601']) == False):
            T = p['TEMPS601'].values
        elif (p['TEMPS602'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPS602']) == False):
            T = p['TEMPS602'].values
        elif (p['TEMPS901'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPS901']) == False):
            T = p['TEMPS901'].values
        elif (p['TEMPS902'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPS902']) == False):
            T = p['TEMPS902'].values
        elif (p['TEMPST01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPST01']) == False):
            T = p['TEMPST01'].values
        elif (p['TEMPRTN1'].isnull().sum() < (0.5*len(p))) & (is_unique(p['TEMPRTN1']) == False):
            T = p['TEMPRTN1'].values
        else:
            #print("Missing temp data: ", castpk)
            stationnodata_T += 1
            continue
        
        dates_np = np.array([datetime.datetime.strptime(x[0:10], '%Y-%m-%d') for x in p['time'].values])
        
        # drop nulls
        T, P, dates_np, nullvals_T = drop_nulls(T, P, dates_np, nullvals_T)
        
        # drop nodata / bad values
        T, P, dates_np, nbadvals_t = drop_bad(-9.99e-29, T, P, dates_np, nbadvals_t)
        T, P, dates_np, nbadvals_t = drop_bad(9.96921E36*-1, T*-1, P, dates_np, nbadvals_t)
        T, P, dates_np, nbadvals_t = drop_bad(9.96921e36*-1, T*-1, P, dates_np, nbadvals_t)
        T, P, dates_np, nbadvals_t = drop_bad(-9.9964, T, P, dates_np, nbadvals_t)
        T, P, dates_np, nbadvals_t = drop_bad(-99.874, T, P, dates_np, nbadvals_t)

        outfile = Rosette_Bottle_Temp(savedir, bottle_str, dates_np, T, P, lat, lon)
        print("Saved", outfile)
        
        n+=1
        n_temp+=1    
        
        #if n > 200:
        #  break
    
    # salinity
    n=0
    stationnodata_P=0
    for idx, row in profiles.iterrows():
        
        castpk = row['profile']
        p = d[d['profile'].values == castpk]
        dep = 0
        
        bottle_str = str(p['profile'].values[0])
        lat = p['latitude'].values[0]
        lon = p['longitude'].values[0]
         
        # depth / pressure
        # sea_water_pressure (dbar), PRESPR01 (decibar), depth (m)
        if (p['sea_water_pressure'].isnull().sum() < (0.5*len(p))) & (is_unique(p['sea_water_pressure']) == False):
            P = p['sea_water_pressure'].values
        elif (p['PRESPR01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['PRESPR01']) == False):
            P = p['PRESPR01'].values
        else:
            #print("Missing pressure data: ", castpk)
            stationnodata_P += 1
            continue
            
        # salinity - check that more than half are populated and not with the same value (indicating placeholder / empty data), convert
        # sea_water_practical_salinity (PSS-78), SSALST01 (PPT), PSALST01 (PSS-78), PSALST02 (PSS-78)
        # FYI PSS_78 is practical salnity / PSU
        if (p['sea_water_practical_salinity'].isnull().sum() == 0) & (is_unique(p['sea_water_practical_salinity']) == False):
            SP = p['sea_water_practical_salinity'].values
        elif (p['PSALST01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['PSALST01']) == False):
            SP = p['PSALST01'].values
        elif (p['PSALST02'].isnull().sum() < (0.5*len(p))) & (is_unique(p['PSALST02']) == False):
            SP = p['PSALST02'].values
        elif (p['SSALST01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['SSALST01']) == False):
            SP = p['SSALST01'].values
            """ Convert PPT to PSU (PSS-78)
            ----------
            https://github.com/TEOS-10/python-gsw/blob/master/gsw/gibbs/practical_salinity.py

            IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
            of seawater - 2010: Calculation and use of thermodynamic properties.
            Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
            UNESCO (English), 196 pp.  See Appendix A.3.
            """
            print("converting from PPT to PSU (PSS-78)")
            SP = (SP - 0.03) * (1.80655 / 1.805)
            SP = np.maximum(SP, 0)  # Ensure that SP is non-negative.
        elif (p['SSALBST01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['SSALBST01']) == False):
            SP = p['SSALBST01'].values
            print("converting from PPT to PSU (PSS-78)")
            SP = (SP - 0.03) * (1.80655 / 1.805)
            SP = np.maximum(SP, 0)  # Ensure that SP is non-negative.
        else:
            #print("Missing salinity: ", castpk)
            stationnodata_S += 1
            continue

        dates_np = np.array([datetime.datetime.strptime(x[0:10], '%Y-%m-%d') for x in p['time'].values])
        
        # drop nulls
        SP, P, dates_np, nullvals_S = drop_nulls(SP, P, dates_np, nullvals_S)

        # drop bad / nodata
        SP, P, dates_np, nbadvals_s = drop_bad(-9.99e-29, SP, P, dates_np, nbadvals_s)
        SP, P, dates_np, nbadvals_s = drop_bad(9.96921E36*-1, SP*-1, P, dates_np, nbadvals_s)
        SP, P, dates_np, nbadvals_s = drop_bad(9.96921e36*-1, SP*-1, P, dates_np, nbadvals_s)
        SP, P, dates_np, nbadvals_s = drop_bad(-9.9964, SP, P, dates_np, nbadvals_s)
        SP, P, dates_np, nbadvals_s = drop_bad(-99.874, SP, P, dates_np, nbadvals_s)

        outfile = Rosette_Bottle_Salin(savedir, bottle_str, dates_np, SP, P, lat, lon)
        print("Saved", outfile)

        n+=1
        n_salin+=1
        
        # REMOVE AFTER DEBUG
        #if n > 200:
        #    break
      
print("Length of file: ", len(d))      
print("Found {} unique rosette bottle grabs.".format(len(profiles)))
print("##########################")
print("Saved {} grabs with SALIN data.".format(n_salin))
print("Skipped {} bad or missing salinity values".format(stationnodata_S))
print("Dropped {} null salinity values".format(nbadvals_s))
print("##########################")
print("Saved {} grabs with TEMP data.".format(n_temp))
print("Skipped {} bad or missing temp values".format(stationnodata_T))
print("Dropped {} null temp values".format(nbadvals_t))
print("##########################")
print("Skipped {} bad or missing pressure values".format(stationnodata_P))

