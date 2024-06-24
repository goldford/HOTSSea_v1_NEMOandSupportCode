
# grabs observations corresponding to a set of model output files
# GO - 20230424
# this script meant to be run on unix (server: graham)
# downloaded 2024-03-05

# 1/ get the list of extracted model data
# 2/ loop through and make a list of corresponding CTD data in OBS-WEST
# 3/ interpolate to model depths
# 4/ export to import glob


import matplotlib.pyplot as plt
import netCDF4
import cmocean as cm
import numpy as np

from scipy import interpolate

import pickle
import os
import sys
import glob
import xarray as xr

sys.path.insert(1, '/project/6006412/goldford/SCRIPTS/')
from utils_GO import new_utils_GO 

shared_p = '/project/6006412/goldford/'
mod_p = 'ANALYSIS/ORAS5-JdF/EXTRACT/CTD/'
mod_out_p = 'ANALYSIS/ORAS5-JdF/CUSTOM-GO/'
obs_p = 'OBS-WEST/CTD/'
obs_out_p = 'OBS-JdF/'
nemo_p = '/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN203/CDF/1990/'


drop_firstlev = True
export_plots = False

# get model depths from random nemo_result (smaller file than mesh_mask)
tmodel = xr.open_dataset(nemo_p + 'SalishSea1500-RUN203_1d_grid_T_y1990m04.nc')
mod_depth = tmodel['deptht']

deplevs = np.zeros(40)
for i in range(0,40):
    deplevs[i] = mod_depth[i]
print(deplevs)

glob_pattern = shared_p + mod_p + "/**/*.pickle"

# Get a list of all model files in the directory and subdirectories
mod_f = glob.glob(glob_pattern, recursive=True)
# sort by date
mod_f.sort(key = lambda x: x[-19:])
time_ar_mod_intrp = np.empty(len(mod_f), dtype='datetime64[s]')
salt_ar_mod_intrp = np.zeros([len(mod_f), 40])
temp_ar_mod_intrp = np.zeros([len(mod_f), 40])

# //////////////////////////////////////////////////////////////
# find obs files to match 
matching_nc_files = []
i=0

for pickle_file in mod_f:
    
    # split the path into components (shared, model specific, subfoldeR)
    subfolder, filename = os.path.split(pickle_file)
    filename_without_ext = os.path.splitext(filename)[0]
    segments = subfolder.split(os.path.sep)

    # Find the indices of the first and last segment corresponding to the input subpaths
    split_shr_p = shared_p.split(os.path.sep)
    sharedp_idx1 = segments.index(split_shr_p[1])
    sharedp_idx2 = segments.index(split_shr_p[len(split_shr_p)-2])
    split_mod_p = mod_p.split(os.path.sep)
    modp_idx1 = segments.index(split_mod_p[0])
    modp_idx2 = segments.index(split_mod_p[len(split_mod_p)-2])

    # path segments for each subpath
    subpath1_segments = segments[1:sharedp_idx2+1]
    subpath2_segments = segments[modp_idx1:modp_idx2+1]
    remaining_segments = segments[modp_idx2+1:]
    
    # Join the subsegments into paths
    subpath1_path = os.path.join(*subpath1_segments)
    subpath2_path = os.path.join(*subpath2_segments)
    remaining_path = os.path.join(*remaining_segments)

    # Construct the expected netCDF file path
    nc_file_path = os.path.join(shared_p, obs_p, remaining_path, filename_without_ext + ".nc")
    if os.path.exists(nc_file_path):
        matching_nc_files.append(nc_file_path)
    
    # //////////////////////////////    
    # interp model data 
    mod_d = pickle.load(open(pickle_file, 'rb'))
    mod_salt = mod_d['salinity']
    mod_temp = mod_d['pTemp']
    mod_time = mod_d['time']
    mod_depth = mod_d['z']
    mod_ndep = len(mod_depth)
    
    t_full, s_full = new_utils_GO.get_model_interpolated_ar(mod_depth, mod_salt, mod_temp, mod_time, deplevs, verbose = False)
    
    salt_ar_mod_intrp[i,:] = s_full
    temp_ar_mod_intrp[i,:] = t_full
    time_ar_mod_intrp[i] = mod_time
    
    if i%200==0:
      print("Interpolating mod file ", i)

    #if i == 5: 
    #  print("T before interp: ", mod_temp)
    #  print("T after interp: ", t_full)
    #  print("S before interp: ", mod_salt)
    #  print("S after interp: ", s_full)

    i+=1

# added by G0 since rarely CTD starts <= 0.5 m
if drop_firstlev == True:
    salt_ar_mod_intrp[:,0] = np.nan
    temp_ar_mod_intrp[:,0] = np.nan

#print("check on interpolated temperature: ", temp_ar_mod_intrp[5,:])
    
print("dumping pickle files for mod data, z-interp.")
pickle.dump(salt_ar_mod_intrp, open(shared_p + mod_out_p + "temp-ORAS5-zinter-semiblind-salt_array.pkl", 'wb'))
pickle.dump(temp_ar_mod_intrp, open(shared_p + mod_out_p + "temp-ORAS5-zinter-semiblind-temp_array.pkl", 'wb'))
pickle.dump(time_ar_mod_intrp, open(shared_p + mod_out_p + "temp-ORAS5-zinter-semiblind-time_array.pkl", 'wb')) 
    

if export_plots == True:
  print("plotting mod data (v interp)")
  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_t = plt.pcolormesh(np.transpose(temp_ar_mod_intrp), vmin = 4, vmax = 15, cmap = cm.cm.thermal)
  plt.colorbar(w_t, label = 'deg C')
  axs.invert_yaxis()
  fig.savefig('temp_JdF_ORAS5_zinterp.jpg')

  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_s = plt.pcolormesh(np.transpose(salt_ar_mod_intrp), vmin = 20, vmax = 33, cmap = cm.cm.thermal)
  plt.colorbar(w_s, label = 'PSU')
  axs.invert_yaxis()
  fig.savefig(shared_p + mod_out_p + 'salt_JdF_ORAS5_zinterp.jpg')

print("Searched for matching obs files for " + str(len(mod_f)) + " mod files")
print("Found this many matching obs files:")
print(len(matching_nc_files))


# with open('./sampleda./ForTereza/prepped_pyapnames/11-03-08_2014h.pickle', 'rb') as pickle_file:
#     model = pickle.load(pickle_file)

# filter out files with all measurements as nan
bad_ctds = []
good_ctds = []
for i in range(0, len(matching_nc_files)):
  try:
    obs = xr.open_dataset(matching_nc_files[i])
    ttime = (obs['time'].values)
    tobs = len(obs['Pres'].values)
    t = np.empty(tobs, dtype='datetime64[s]')
    t[:] = ttime[0]
    if bool(obs['cTemp'].isnull().all()) is not True:
        good_ctds.append(matching_nc_files[i])
    else:
        print("found cast w/ all nans for temp")
  except:
    bad_ctds.append(matching_nc_files[i])
  
print("tried to open this many files: ", i)
print("didn't work for this many files: ", len(bad_ctds))

if export_plots == True:
  print("plotting obs data.")
  # check for empty records and plot obs data as hovmoller
  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')

  for i in range(0, len(good_ctds)):
    obs = xr.open_dataset(good_ctds[i])
    ttime = (obs['time'].values)
    tobs = len(obs['Pres'].values)
    t = np.empty(tobs, dtype='datetime64[s]')
    t[:] = ttime[0]
    w=axs.scatter(t, obs['Pres'].values, c = obs['cTemp'].values, s = 2, vmin = 4, vmax = 15, cmap = cm.cm.thermal)

  plt.colorbar(w, label = 'deg C')
  axs.invert_yaxis()

  fig.savefig(shared_p + mod_out_p + 'raw_ctd_JdF.jpg')

# sort them by date
ctd_bydate = good_ctds.copy()
# sort based on the date in the final positions of file name
ctd_bydate.sort(key = lambda x: x[-19:])
#print("first ctds to check sorting worked: ", ctd_bydate[0:10])


time_array = np.empty(len(ctd_bydate), dtype='datetime64[s]')
salt_array = np.zeros([len(ctd_bydate), 40])
temp_array = np.zeros([len(ctd_bydate), 40])

for i in range(0,len(ctd_bydate)):

  tstr = ctd_bydate[i]
  t_full, s_full, ttime  = new_utils_GO.get_model_interpolated(tstr, deplevs, verbose = False)

  salt_array[i,:] = s_full
  temp_array[i,:] = t_full
  time_array[i] = ttime
    
  if i%200==0:
    print("Interpolating obs file ", i)

# added by G0 since rarely CTD starts <= 0.5 m
if drop_firstlev == True:
  salt_array[:,0] = np.nan
  temp_array[:,0] = np.nan           

if export_plots == True:
  
  print("plotting obs data (v interp)...")
  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_t = plt.pcolormesh(np.transpose(temp_array), vmin = 4, vmax = 15, cmap = cm.cm.thermal)
  plt.colorbar(w_t, label = 'deg C')
  axs.invert_yaxis()
  fig.savefig(shared_p + mod_out_p +'temp_JdF_OBS_zinterp.jpg')

  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_s = plt.pcolormesh(np.transpose(salt_array), vmin = 20, vmax = 33, cmap = cm.cm.thermal)
  plt.colorbar(w_s, label = 'PSU')
  axs.invert_yaxis()
  fig.savefig(shared_p + mod_out_p +'salt_JdF_OBS_zinterp.jpg')
  

print("dumping pickle files for obs data, z-interp.")
pickle.dump(salt_array, open(shared_p + obs_out_p + "CTD-salt_array.pkl", 'wb'))
pickle.dump(temp_array, open(shared_p + obs_out_p + "CTD-temp_array.pkl", 'wb'))
pickle.dump(time_array, open(shared_p + obs_out_p + "CTD-time_array.pkl", 'wb')) 

print("finished")
