''' convert DFO CTD data
    G Oldford Jul 2021

    Based on convert_HakaiCTDs.py
    data in: large CSV of DFO IOS CTD casts
    data out: many .NC files (one per cast)

    Source data from 1978 to 2020 downloaded from CIOOS Pacific
    https://data.cioospacific.ca/erddap/tabledap/IOS_CTD_Profiles.html
    Used Salish Sea bounding box but this included stations off WCVI
'''

import glob
import os
import gsw
import numpy as np
import pandas as pd
import analysispkg
import re

import sys
sys.path.insert(1, '/project/6006412/goldford/py3-8-10_Jul9/lib/python3.8/site-packages/python-analysis-package/etc/data_conversion/')
from data_conversion_functions import write_to_netcdf

indir = '/project/6006412/goldford/data/evaluation/IOS_CTD_Profiles/'
savedir = '/project/6006412/goldford/data/evaluation/IOS_CTD_Profiles/prepped/'

os.makedirs(savedir, exist_ok=True)

files = sorted(glob.glob(os.path.join(indir, "*.csv")))
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
    d = pd.read_csv(file, skiprows=range(1,2), dtype={"instrument_model": object, 
                                                      "instrument_serial_number": object, 
                                                      "latitude": float, "longitude": float,
                                                      "depth": float, "CNDCST01": float,
                                                      "DOXMZZ01": float, "DOXYZZ01": float, 
                                                      "PRESPR01": float, "PSALST01": float,
                                                      "PSALST02": float, "sea_water_practical_salinity": float, 
                                                      "sea_water_pressure": float, "sea_water_temperature": float,
                                                      "SSALST01": float, "TEMPS601": float, "TEMPS602": float, 
                                                      "TEMPS901": float, "TEMPS902": float, "TEMPST01": float
                                                      })
    print(d.columns)
    profiles = d[['profile']].drop_duplicates()
    print("Unique casts: ", len(profiles))
    print("Length of file: ", len(d))
    print(d.head())
   
    # first file row contains vars, units in second row
    d_units = pd.read_csv(file, nrows=1)
    print("Units: ", d_units)
    
    for idx, row in profiles.iterrows():
        
        castpk = row['profile']
        p = d[d['profile'].values == castpk]
        dep = 0
        
        #stationstr = p['geographic_area'].values[0]
        stationstr = str(p['project'].values[0])
        datetimestr = p['time'].values[0]
        lat = p['latitude'].values[0]
        lon = p['longitude'].values[0]
         
        # temperature
        # metadata missing info but indicates ITS90 temp thus in situ not potential.
        # checks: not all values the same, not > 50% values are nulls with priority given
        #         to sea_water_temperature field which appears to be a standardized field 
        # https://data.cioospacific.ca/erddap/info/IOS_CTD_Profiles/index.html
        # -99.847 is used as null value in some records
        # sea_water_temperature (degC), TEMPS601 (degC), TEMPS602 (degC), TEMPS901 (degC), TEMPS902 (degC), TEMPST01 (degC)
        if (p['sea_water_temperature'].isnull().sum() > (0.5*len(p))) & (is_unique(p['sea_water_temperature']) == False):
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
        else:
            print("No complete temp data found in cast: ", castpk)
            #T = np.empty([len(p)])
            continue
        
        # salinity
        # sea_water_practical_salinity (PSS-78), SSALST01 (PPT), PSALST01 (PSS-78), PSALST02 (PSS-78)
        # PSS_78 is practical salnity / PSU
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
            SP = (SP - 0.03) * (1.80655 / 1.805)
            SP = np.maximum(SP, 0)  # Ensure that SP is non-negative.
        else:
            print("No complete salinity data found in cast: ", castpk)
            continue
            #SP = np.empty([len(p)])

        # depth / pressure
        # sea_water_pressure (dbar), PRESPR01 (decibar), depth (m)
        if (p['sea_water_pressure'].isnull().sum() < (0.5*len(p))) & (is_unique(p['sea_water_pressure']) == False):
            P = p['sea_water_pressure'].values
        elif (p['PRESPR01'].isnull().sum() < (0.5*len(p))) & (is_unique(p['PRESPR01']) == False):
            P = p['PRESPR01'].values
        else:
            print("No complete pressure data found in cast: ", castpk)
            #P = np.empty([len(p)])
            continue

        # -9.99E-29 appears to be a flag value for missing/bad data in Hakai CTD data
        T,SP,P,badvals = replace_bad(-9.99e-29, T,SP,P,badvals)
        # bad value for first layer of salinity in DFO data
        T,SP,P,badvals = replace_bad(-9.9964, T,SP,P,badvals)
        T,SP,P,badvals = replace_bad(-99.874, T,SP,P,badvals)
        # to do: flag and remove nulls as bad vals?

        absSal = gsw.conversions.SA_from_SP(SP,P,lon,lat) # abs salin
        consTemp = gsw.conversions.CT_from_t(absSal, T, P) # conservative temp
        pDen = gsw.rho(absSal, consTemp, P) # dens 
        
        # replace special chars
        replaced_stationstr = re.sub('[,. &]', '_', stationstr)
        
        outfile = write_to_netcdf().CTD_Cast(savedir, replaced_stationstr, datetimestr, P, T, consTemp, SP, absSal, pDen, dep, lat, lon)
        print("Saved", outfile)

        n+=1
        
        # REMOVE AFTER DEBUG
        #if n > 200:
        #    break

print("Saved {} profiles.".format(n))
print("Dropped {} bad values (-9.99E-29 or -9.9964, -99.874)".format(badvals))


