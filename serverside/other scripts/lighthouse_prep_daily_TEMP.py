
''' convert DFO Lighthouse Daily data
    G Oldford
    Last modified: Dec 11 2021
    
    left off w/ line 254, salinity temp data
    
    Based on code from FA12 project convert_Hakai_CTDs
    data in: CSV of DFO Lighthouse Daily
    data out: many .NC files (one per lightstation and separated by temp, salin)
    
    To do: 
    - 1) get lighthouse timestamp from TG data, parse, replace
    - 2) write salinity and temp to NC file, one per lighthouse (don't separate temp and salinity this time)
    - 3) the closest TG link table needs to be expanded given patchy data. particularly TG gauges 7917, 7982
    - 4) 7532 seems to be mislabelled in the link file - either the TG is Sand Heads or it's Whaler Bay (7532) that's closest to active pass

    Source data from 1978 to 2020 downloaded from CIOOS Pacific
    https://data.cioospacific.ca/erddap/tabledap/BCSOP_daily.html
    metadata:
    https://catalogue.cioos.ca/dataset/ca-cioos_654a4ece-7271-4f5a-ba60-b50a62dbd051
    https://data.cioospacific.ca/erddap/info/BCSOP_daily/index.html
     
    
    # pseudocode
    # for each record in the lighthouse CSV file
    #   get the ID of the nearest tidal gauge station with data 
    #   get the TG data itself 
    #   find the time between 6 am and 6 pm w/ highest tide
    #   convert the time 
    #   overwrite the lighthouse time stamp
'''
import netCDF4 as nc
import glob
import os
from os.path import abspath, expanduser, exists, join
import gsw
import numpy as np
import pandas as pd
import datetime as datetime
import time
import csv
import re
import math
#import analysispkg

indir_lh = '/project/6006412/goldford/data/evaluation/Lighthouses/'
indir_tg = '/project/6006412/mdunphy/TideGauges/'
savedir = '/project/6006412/goldford/data/evaluation/Lighthouses/processed/'

stationlist = 'CHS_PAC_stationlisting.txt'  # pipe | separated
#STATION_NUMBER|STATION_NAME|LATITUDE|A.LONGITUDE*(-1)|MIN(D.YEAR)|MAX(D.YEAR)
#7020|Sooke , B.C.|48.3695|-123.726|1972|2011

tg_template = '{}_HOURLY.DAT' # separated by spaces (inconsistent)
# 1972/12/01 00:00 999.999
# 1972/12/01 01:00   1.990

os.makedirs(savedir, exist_ok=True)

#files_tg = sorted(glob.glob(os.path.join(indir_tg, "*.DAT")))
file_lh_tg = indir_lh + "LighthousesAndClosestTG.csv"


# assign each lighthouse a numeric ID 
def get_LH_ID (LH_ID_dict):

  # IDs match those in LighthousesAndClosesTG
  my_dict = dict([(0,'active_pass_lightstation'), 
                  (1,'amphitrite_point_lightstation'),
                  (2, 'cape_beale_lightstation'),
                  (3, 'cape_mudge_lightstation'),
                  (4, 'chrome_island_lightstation'),
                  (5, 'departure_bay_(pbs)'),
                  (6, 'entrance_island_lightstation'),
                  (7, 'race_rocks_lightstation'),
                  (8, 'sheringham_point_lightstation'),
                  (9, 'sisters_islets_lightstation'),
                  (10, 'west_van_labs_(caer)'),
                  (11, 'bamfield_marine_sciences_centre')
                  ])
  #to do: finish this. rather do a loop over dict than have it here -GO
  LH_ID = 0
  return LH_ID

# finds closest TG with data for given date and attempts to retrieve data
def get_TG_time (LH_key, df_LH_TG, LH_date):
# link table
  # FID,lat,lon,LH_ID,name,ClosestTG,TG_ID,MIN_YEAR,MAX_YEAR,rank_close
  # 0,48.2979,-123.5316,7,race_rocks_lightstation,"Pedder Bay, B.C.",7080,1967,1969,1
  
  LH_year = LH_date.year
  
  # get the TG's for lighthouse (LH_key) for year 
  df_LH_TG_filtered = df_LH_TG[df_LH_TG['LH_ID'].values == LH_key]
  df_LH_TG_filtered = df_LH_TG_filtered[(df_LH_TG_filtered['MIN_YEAR'].values <= LH_year) & (df_LH_TG_filtered['MAX_YEAR'].values >= LH_year)]
  
  if len(df_LH_TG_filtered) > 0:
    df_LH_TG_filtered.sort_values(by=['rank_close'])

    for index, row in df_LH_TG_filtered.iterrows():
      #print(row['TG_ID'])
      #print(row['rank_close'])
      TG_ID = row['TG_ID']
      
      TG_hightide_time = get_TG_data(TG_ID, LH_date)

      if TG_hightide_time != 999:
        break

  else:
    print("No matching tidal gauge found for LH " + str(LH_key))
    TG_found = False
    TG_ID = 0
    TG_hightide_time = 999
    
  return TG_hightide_time

# TG data stored in one file per station (09476_HOURLY.DAT)
# TG DAT file has header 24 lines long, no column headings, nodata= 999.999
# 2018/12/05 06:00 999.999
# 2018/12/05 18:00   4.586 note: three spaces (9 digits for heights)
def get_TG_data (TG_ID, LH_date):

  # try to open dat file
  TG_dat_df = pd.read_csv(indir_tg + str(0) + str(TG_ID) + "_HOURLY.DAT", delim_whitespace=True, skiprows=24)
  
  # produce two cols with python dates
  TG_dat_df['date_day'] = TG_dat_df.iloc[:,0].astype(str)
  TG_dat_df['date_day_py'] = pd.to_datetime(TG_dat_df['date_day'])
  TG_dat_df['date_day_py'] = pd.to_datetime(TG_dat_df['date_day'])
  
  TG_dat_df['hour'] = TG_dat_df.iloc[:,1].astype(str)
  TG_dat_df['date_day_hr'] = TG_dat_df.loc[:,'date_day'] + " " + TG_dat_df.loc[:,'hour']
  TG_dat_df['date_day_hr_py'] = pd.to_datetime(TG_dat_df['date_day_hr'])
  
  TG_dat_df['ssh'] = TG_dat_df.iloc[:,2].astype(float)
  
  TG_dat_df_filtered = TG_dat_df[(TG_dat_df['date_day_py'].dt.year == LH_date.year) & 
                                 (TG_dat_df['date_day_py'].dt.month == LH_date.month) & 
                                 (TG_dat_df['date_day_py'].dt.day == LH_date.day)]
  
  # loop through the 0700 to 1700 and get the high tide
  ssh = 0
  hightide_time = 0
  if len(TG_dat_df_filtered) > 0: 
    for index, row in TG_dat_df_filtered.iterrows():
      if (row['date_day_hr_py'].hour > 6) & (row['date_day_hr_py'].hour < 17):
        if (row['ssh'] > ssh) & (row['ssh'] < 999): 
          ssh = row['ssh']
          hightide_time = row['date_day_hr_py']
          
    if ssh == 0: 
      hightide_time = 999
  else: 
    hightide_time = 999

  return hightide_time

#unused - just template code
def update_NC_timestamp (timevar, LHfile):
    
     # sample code from MD adjustERA5 script  --->
    with nc.Dataset(infile1) as ncf:
        tmp = np.zeros(ncf[var].shape) - 999
        tmp[:-1,:,:] = ncf[var][1:,:,:]  # second record becomes first, etc
    if year == year_end:
        # we don't have a 2021 file, duplicate last record as a workaround
        tmp[-1, :, :] = tmp[-2, :, :]
    else:
        with nc.Dataset(infile2) as ncf:
            tmp[-1,:,:] = ncf[var][0,:,:]  # first record from next year is last record of this year
    return tmp

# added feb 11 - GO
def makePath(path):
    return abspath(expanduser(path))

# added feb 11 - GO
def BC_Lighthouse_Daily(saveDir, station, dates_, Ts_, Ss_, lat_, lon_, depths_,yr):
    writeName = station.replace('/', '-') + '_SalinTemp_' + str(yr) + '.nc'
    
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
        time_counter[:] = nc.date2num(dates_, time_counter.units)

        Temp = f.createVariable('temperature', 'f8', 'time_counter', fill_value=-999)
        Temp[:] = Ts_
        Temp.type = 'in-situ'
        Temp.units = 'degC'
        Temp.description = 'recorded temperature'
        
        Sal = f.createVariable('salinity', 'f8', 'time_counter', fill_value=-999)
        Sal[:] = Ss_
        Sal.type = 'practical'
        Sal.units = 'psu'
        Sal.description = 'recorded salinity'
        
       # cT = f.createVariable('cTemp', 'f8', 'time_counter')
       # cT[:] = consTemp
       # cT.units = 'degC'
       # cT.description = 'conservative temperature'

        Lat = f.createVariable('latitude', 'f4', 'single')
        Lat[:] = lat_
        Lat.units = 'degrees north'
        Lat.description = 'cast latitude'

        Lon = f.createVariable('longitude', 'f4', 'single')
        Lon[:] = lon_
        Lon.units = 'degrees east'
        Lon.description = 'cast longitude (-180 -> 180 degrees)'
        
        Dep = f.createVariable('depth', 'f4', 'time_counter')
        Dep[:] = depths_
        Dep.units = 'metres'
        Dep.description = 'approximate depth in metres of measurement'

        return svPath

# =====================================================
# main loop

# IDs match those in LighthousesAndClosestTG
my_dict = dict([#(0,'active_pass_lightstation'), 
                #(1,'amphitrite_point_lightstation'),
                #(2, 'cape_beale_lightstation'),
                #(3, 'cape_mudge_lightstation'),
                #(4, 'chrome_island_lightstation'),
                #(5, 'departure_bay_(pbs)'),
                #(6, 'entrance_island_lightstation'),
                #(7, 'race_rocks_lightstation'),
                #(8, 'sheringham_point_lightstation'),
                (9, 'sisters_islets_lightstation')
                #(10, 'west_van_labs_(caer)'),
                #(11, 'bamfield_marine_sciences_centre')
                  ])

# link table
# FID,lat,lon,LH_ID,name,ClosestTG,TG_ID,MIN_YEAR,MAX_YEAR,rank_close
# 0,48.2979,-123.5316,7,race_rocks_lightstation,"Pedder Bay, B.C.",7080,1967,1969,1
df_lh_tg = pd.read_csv (indir_lh + "LighthousesAndClosestTG.csv")

# lighthouse data
LHfile = indir_lh + "original/BCSOP_daily_SalishSea_1978to2019.csv"
lh_df = pd.read_csv(LHfile, skiprows=range(1,2), 
                      dtype={"profile": object,
                      "latitude": float, "longitude": float,
                      "TEMPTC01": float, "PSALPR01": float,},
                      usecols=['profile','latitude','longitude','TEMPTC01','PSALPR01','time'], 
                      low_memory=False
                      )
                      
print("Reading Lighthouse data, length: " + str(len(lh_df)))
lh_df['date'] = pd.to_datetime(lh_df['time'])
lh_df['year'] = lh_df['date'].dt.year
                        
# for each lighthouse in dictionary get corresponding data record
# get date and parse to pydate
# get closest Tidal Gauge with data for that date
# retrieve TG data and determine time of day with high tide
# return high tide as pydate
# get salinity, temp values
# write to file

success = 0
fail = 0

for key, lighthouse in my_dict.items():

  lighthousestr = str(lighthouse)
  print("retrieving tidal gauge data closest to lighthouse " + str(lighthouse))
  
  # if I don't do by year the process gets killed
  for yr in range(2006,2007):
  
   # data
   lh_df_one = lh_df[(lh_df['profile'].values == lighthouse) & (lh_df['year'].values == yr)]
   
   if len(lh_df_one)==0:
     print("no records found")
     continue
   
   lat = lh_df_one['latitude'].values[0]
   lon = lh_df_one['longitude'].values[0]
   
   i = 0 
   Ts = []
   Ss = []
   depths = []
   TG_hightide_times = []
 
   for index, row in lh_df_one.iterrows():
     
     lh_date = row['date']
 
     # debug
     #if lh_date.year != 2006:
     #  continue
     
     # profile = lighthouse name
     #print(row['profile'], row['date'])
     #print(row)
     
     # get new date / time of day from TG
     TG_hightide_time = get_TG_time(key, df_lh_tg, lh_date)
     if TG_hightide_time == 999:
       fail += 1
       #print("no high tide time")
       continue
     else:
       #print("high tide time")
       success += 1
       TG_hightide_times.append(TG_hightide_time)
     
     print(TG_hightide_time)
 
     # 2022-02-11
     T = row['TEMPTC01']
     # error catching
     try:
       if T is None: 
         print('T is None')
         T = -999
     except NameError:
         print ("T is not defined")
         T = -999
     else:
         if math.isnan(T):
           #print ("T is nan")
           T = -999     
     Ts.append(T)
     
     # salinity
     S = row['PSALPR01']
     try:
       if S is None: 
         print('S is None')
         S = -999
     except NameError:
         print ("S is not defined")
         S = -999
     else:
         if math.isnan(S):
           #print ("S is nan")
           S = -999
     #print(S)
     Ss.append(S)
     
     # assumed
     dep=0.5
     depths.append(dep)
 
     # debug exitor
     #i += 1
     #if i > 2:
     #  break   
 
   if len(Ts) > 0:
     outfile = BC_Lighthouse_Daily(savedir, lighthousestr, TG_hightide_times, Ts, Ss, lat, lon, depths, yr)
   else: 
     "no records found, continuing"
     continue
   print("Saved", outfile)


print(str(success/(success+fail)*100) + "% success at finding high tide time for lighthouse data from tide gauges")
  
  
  
  # scrap
  #with nc.Dataset(LHfile) as ncf:
  #  tmp = np.zeros(ncf[var].shape) - 999
  #  tmp[:-1,:,:] = ncf[var][1:,:,:]  # second record becomes first, etc
  
  
  
  
  