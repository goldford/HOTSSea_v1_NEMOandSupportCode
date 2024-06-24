''' convert Nanoose Stn (military) files 2002 -on
    G Oldford Aug 2022
    
    Based on code in python-analysis-package etc/data_conversion_functions.py 
    data in: many CTD files in many folders from Nanoose Station 
    data out: many .NC files (one per CTD cast measurement)

'''
import glob
import os
import gsw
import numpy as np
import pandas as pd
#import analysispkg
import re
import netCDF4 as nc
import datetime as datetime
import sys

# utils script
sys.path.insert(1, '/project/6006412/goldford/data/evaluation/')
from utils_GO import write_to_netcdf_GO 
from utils_GO import general_utils

#sys.path.insert(1, '/project/6006412/goldford/pypkgtemp/python-analysis-package/etc/data_conversion/')
#from data_conversion_functions import write_to_netcdf

#indir = '/project/6006412/goldford/data/evaluation/nanoose_stn/wprop_2002_2018_ctd/**/'
indir = '/project/6006412/goldford/data/evaluation/nanoose_stn/wprop_2002_2018_ctd/'
#indir = '/project/6006412/goldford/data/evaluation/nanoose_stn/NANOOSE/9679/'
saveDir = '/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_2003on/'
os.makedirs(saveDir, exist_ok=True)

files_MED = sorted(glob.glob(os.path.join(indir, "*.MED"), recursive=True))
#files_CTD = sorted(glob.glob(os.path.join(indir, "*.CTD"), recursive=True))
files_CTD = sorted(glob.glob(os.path.join(indir, "*.ctd"), recursive=True))

print("found MED files: ", len(files_MED))
print("found CUR files: ", len(files_CTD))

n=0
badvals=0

fi = 0
cnt = 0 

## main loop 
for file in files_CTD: 
    # for debug      
    fi = fi + 1
    #if fi > 10:
    #  break
    
    print("Working on ", file)
    
    with open(file,'r',encoding='latin-1') as f:
      lines = [line[:-1] for line in f] # note may need line[:-1] or line[:-2] to truncate /n -GO
    
    timeRange = []    # my code - GO
    # flags for metadata about vars
    st_check = 0; et_check = 0; ti_check = 0
    sn_check = 0; wd_check = 0; lt_check = 0
    ln_check = 0; id_check = 0; eh_check = 0
    old_file = -1; mn_check = 0; mission = ""
    
    nameLines = [] # sample code begin -GO
    for i, line in enumerate(lines): # loop of header info
        if 'EVENT NUMBER' in line:

            if ':' in line:
              l_i = line.find(':') + 1
              station1 = line[l_i:].replace(' ', '')
              sn_check = 1
              old_file = 1
              #print("Station (event num):", station1)
        
        if 'MISSION' in line:

            if ':' in line:
              l_i = line.find(':') + 1
              mission = line[l_i:].replace(' ', '')
              mn_check = 1
              #print("Station (event num):", station1)
              
        if 'TIME INCREMENT' in line:
            subLine = line[line.find(':') + 2:line.rfind('!')]
            floats = np.asarray(subLine.split()).astype(float)
            dt = floats[0]*86400. + floats[1]*3600. + floats[2]*60. + floats[3] + floats[4]/1000.
            
            ti_check = 1
            old_file = 0
            #print("got time increment.")

        if len(timeRange) < 2:
            if 'START TIME' in line:
            
              #print("found START TIME.")
            
              #sT = line[line.find('UTC') + 4:line.rfind('.')]
              tz=0
              if line.find('GMT')!=-1:
                tz = line.find('GMT')
              elif line.find('gmt')!=-1:
                tz = line.find('gmt')
              elif line.find('PST')!=-1:
                tz = line.find('PST')              
              elif line.find('PDT')!=-1:
                tz = line.find('PDT')
              elif line.find('pst')!=-1:
                tz = line.find('pst')
              elif line.find('UTC')!=-1:
                tz = line.find('UTC')
              elif line.find('utc')!=-1:
                tz = line.find('utc')
              else:
                print("no time zone found for start time, skipping.")
                break
              
              sT = line[tz + 4:line.rfind('.')]
              sT = datetime.datetime.strptime(sT,'%Y/%m/%d %H:%M:%S')
              timeRange.append(sT)
              
              st_check = 1
              print(sT)

        if 'LATITUDE' in line:
            subLine = line[line.find(':') + 1:line.rfind('!')]
            lat = np.asarray(subLine.split()[:2]).astype(float)
            lat = lat[0] + lat[1]/60.
            
            lt_check = 1
            #print("converted latitude.", lat)

        if 'LONGITUDE' in line:
            subLine = line[line.find(':') + 1:line.rfind('!')]
            lon = np.asarray(subLine.split()[:2]).astype(float)
            lon = lon[0] + lon[1]/60.

            if 'W' in line:
              s_subLine = subLine.split()
              if len(s_subLine) > 2:
                if s_subLine[2] == 'W':
                  lon = lon * -1.
              else:
                if s_subLine[1] == 'W':
                  lon = lon * -1.
            
            ln_check = 1
            #print("converted longitude: ", lon)
            
        if '*END OF HEADER' in line:
            skipheader = i + 1
            eh_check = 1
            break

        # check if different ctd file format
        if 'INITIAL_LATITUDE' in line:
            lat1 = float(line[line.find('=') + 1:line.find(',')])
            print("warning different lat format detected: ", lat1)
        if 'INITIAL_LONGITUDE' in line:
            lon1 = float(line[line.find('=') + 1:line.find(',')])
            print("warning different lon format detected: ", lon1)
        if 'END_LATITUDE' in line:
            lat2 = float(line[line.find('=') + 1:line.find(',')])
            print("warning different lat format detected: ", lat2) 
        if 'END_LONGITUDE' in line:
            lon2 = float(line[line.find('=') + 1:line.find(',')])
            pprint("warning different lon format detected: ", lon2)
        if '-- DATA --' in line:
            skipHeader = i + 1
            print("warning alternative format found. see data line.")
            break
        if 'MAX_DEPTH' in line:
            dep = float(line[line.find('=') + 1:line.find(',')])
            print("warning alt dep format found:", dep)
        if 'START_DATE_TIME' in line:
            print("warning alternative format found. see start_date_time.")
        if 'EVENT_QUALIFIER1' in line:
            station2 = line[line.find('=') + 1:line.find(',')]
            print("warning. event qualifier found: ", station2)
            


    # //////////////////////////////////////
    # cross check minimum info is found
    
    # exception
#    if wd_check == 0:
#      print("no water dep found. continuing anyway.")
#      dep = 0; wd_check = 1
#    if id_check == 0:
#      print("no inst dep found.  continuing anyway.")   
#    if et_check == 0:
#      print("no end time found.continuing anyway.")   
#    if ti_check == 0:
#      print("no time inc found. continuing anyway.")  
#    if wd_check == 0:
#      print("no water dep found. continuing anyway.") 
      
    if st_check == 0:
      print("no start time found. skipping.")
      continue 
    elif sn_check == 0:
      if mn_check == 1:
        station1 = mission
      else:
        print("no stn name found. skipping.")
        continue
    elif lt_check == 0:
      print("no lat found. skipping.")
      continue
    elif ln_check == 0:
      print("no lon found. skipping.")
      continue
    elif eh_check == 0:
      print("no end of header found.  skipping.")
      continue
    elif old_file == 0:
      print("time increment field found incidating newer file format. skipping for now.")
      continue
    else:
      print("basic metadata found!")
    
#    if lat2 < 0.:
#        lat2 = lat1
#    if lon2 < -180.:
#        lon2 = lon1
#    lat = 0.5 * (lat1 + lat2)
#    lon = 0.5 * (lon1 + lon2)
    
    #print("lat:", lat)
    #print("lon:", lon)
#
    
    # ////////////////////////////////////////
    # get index + confirm existence of channels
    TChan = -1; SChan = -1; PChan = -1;
    # some stations/insts have a 'low' and 'high' temp
    TChan_low = -1; TChan_high = -1;
    
    startRead = 0
    for lN,line in enumerate(lines):
      if startRead:
        
        if TChan == -1:
        
          if ('Temperature' in line) and (':Low_Res' not in line) and (':High_Res' not in line):
            TChan = int(line.split()[0]) - 1
            #print("found temp:", TChan)
          elif ('Temperature' in line) and (':Low_Res' in line):
            TChan_low = int(line.split()[0]) - 1
            #print("found temp_low")
          elif ('Temperature' in line) and (':High_Res' in line):
            TChan_high = int(line.split()[0]) - 1
            #print("found temp_high")
            
        if SChan == -1:
          if 'Salinity' in line:
            SChan = int(line.split()[0]) - 1
            #print("found salinity:", SChan)
            
        if PChan == -1: 
          if 'Pressure' in line:
            PChan = int(line.split()[0]) - 1
            #print("found pressure: ", PChan)
        
      if '$TABLE: CHANNELS' in line:
        startRead = 1
      if '$END' in line:
        break

    if (TChan == -1) and (TChan_low == -1) and (TChan_high == -1):
      print("no temp channel. skipping.")
      continue 
    elif (SChan == -1):
      print("no salin channel. skipping.")
      continue 
    elif (PChan == -1):
      print("no pressure channel. skipping.")
      continue 
    else:
      print("all channels found!")
    
    
    # /////////////////////
    # Get data
    T,S,P = np.array([]), np.array([]), np.array([])
    
    lmk = 0
    
    # process all records, remove spaces
    for n,line in enumerate(lines[skipheader:]):
      entries = []
      for e in line.split(' '):
        
        if e != '':
          if '-' in e:
            l = e.split('-')
            for i,l1 in enumerate(l):
              if i == 0:
                entries.append(l1)
              else:
                entries.append('-' + l1)
          else:
             # sometimes overflow of digits causes issues
             # - assumes max width col is 7 digits
             # e.g., '9.367287.8'  
            if (len(e) > 7) and (str(e).count('.') > 1):
              l1 = e[:6]; l2 = e[6:]
              #print("l1: ", l1)
              #print("l2: ", l2)
              entries.append(l1)
              entries.append(l2)
            else:
              entries.append(e)

      for e in entries:
        if e == '':
          entries.remove(e)
        if e == ' ':
          entries.remove(e)

      if len(entries) < 2: 
        continue
        
      if PChan != -1: 
        P = np.hstack([P,float(entries[PChan])])
      else:
        print("no P found. skipping.")
        continue

      if (TChan != -1):
        T = np.hstack([T,float(entries[TChan])])
      elif (TChan == -1) and (TChan_low == -1) and (TChan_high != -1):
        T = np.hstack([T,float(entries[TChan_high])])
      elif (TChan == -1) and (TChan_low != -1) and (TChan_high == -1):
        T = np.hstack([T,float(entries[TChan_low])])
      elif (TChan == -1) and (TChan_low != -1) and (TChan_high != -1):
        T = (np.hstack([T,float(entries[TChan_high])]) + np.hstack([T,float(entries[TChan_low])]))/ 2 
      else:
        print("no T found. skipping.")
        continue

      if SChan != -1: 
        S = np.hstack([S,float(entries[SChan])])
      else:
        print("no S found. skipping.")
        continue
  
    # replace -9999 
    flags = {'P': 0, 'T': 0, 'S': 0}
    
    # invalidate if
    P[P < 0.] = np.nan
    P[P > 4000.] = np.nan
    T[T < -2.] = np.nan
    T[T > 40.] = np.nan
    S[S < 0.] = np.nan
    S[S > 60.] = np.nan

    if all(np.isnan(S)):
        flags['S'] = False
        print("warning. S's are all nans. skipping.")
        continue
    if all(np.isnan(T)):
        flags['T'] = False  
        print("warning.T's are all nans. skipping.")
        continue  
    if all(np.isnan(P)):
        flags['P'] = False
        print("warning. P's are all nans. skipping.")
        continue
    
    T, S, P, badvals = general_utils.remove_nan(T, S, P, badvals)
#      
#    print("found P's - len: ", len(P))
#    print("found T's - len: ", len(T))
#    print("found S's - len: ", len(S))
    
    absSal = gsw.conversions.SA_from_SP(S,P,lon,lat) # abs salin
    consTemp = gsw.conversions.CT_from_t(absSal, T, P) # conservative temp
    pDen = gsw.rho(absSal, consTemp, P) # dens 
    
#    print("found absSal's - len: ", len(absSal))
#    print("found consTemp: ", len(consTemp))
#    print("found pDen: ", len(pDen))
#    print("mean temp: ", np.mean(T))
#    print("mean Salin: ", np.mean(S))
    
    
    dep = 0 # water depth at point of cast
    
    stationstr = re.sub('[,. &]', '_', station1)
    stationstr = "Nanoose_" + stationstr
    
    outfile = write_to_netcdf_GO().nanoose_CTD_Cast(saveDir, stationstr, sT, P, T, consTemp, S, absSal, pDen, dep, lat, lon)
    print("success converting", outfile)
    cnt=cnt+1
    
print("total files attempted:", fi)
print("total files succeeded:", cnt)
   