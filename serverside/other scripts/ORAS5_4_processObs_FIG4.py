
# compares pre-extracted pre-interpolated files from CTDs near mouth of JdF
# to the ORAS5 ocean model files for same area

import pandas as pd
import datetime

import matplotlib.pyplot as plt
import netCDF4
import cmocean as cm
import numpy as np

import pickle
import os
import glob
import xarray as xr

import warnings

shared_p = '/project/6006412/goldford/'
mod_p = 'ANALYSIS/ORAS5-JdF/EXTRACT/'
obs_p = 'OBS-JdF/'
nemo_p = '/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN203/CDF/1990/'
out_p = 'ANALYSIS/ORAS5-JdF/CUSTOM-GO/'

mod_f_s = 'temp-ORAS5-zinter-semiblind-salt_array.pkl'
mod_f_t = 'temp-ORAS5-zinter-semiblind-temp_array.pkl'
mod_f_ti = 'temp-ORAS5-zinter-semiblind-time_array.pkl'
obs_f_s = 'CTD-salt_array.pkl'
obs_f_t = 'CTD-temp_array.pkl'
obs_f_ti = 'CTD-time_array.pkl'

# open the delicious pickles 
# obs data interpolated
obs_temp = pickle.load(open(shared_p + obs_p + obs_f_t, 'rb'))
obs_salt = pickle.load(open(shared_p + obs_p + obs_f_s, 'rb'))
obs_time = pickle.load(open(shared_p + obs_p + obs_f_ti, 'rb'))

# mod data interpolated
mod_temp = pickle.load(open(shared_p + out_p + mod_f_t, 'rb'))
mod_salt = pickle.load(open(shared_p + out_p + mod_f_s, 'rb'))
mod_time = pickle.load(open(shared_p + out_p + mod_f_ti, 'rb'))

# get deplevs for plotting
tmodel = xr.open_dataset(nemo_p + 'SalishSea1500-RUN203_1d_grid_T_y1990m04.nc')
mod_depth = tmodel['deptht']
deplevs = np.zeros(40)
for i in range(0,40):
    deplevs[i] = mod_depth[i]

# shapes don't match between mod / obs yet!
print("Shapes below probably don't match due to some CTDs with all nans dropped --> ")
print(mod_time.shape)
print(obs_time.shape)
print(mod_salt.shape)
print(mod_temp.shape)
print(obs_salt.shape)
print(obs_temp.shape)

# Extract year and day from datetime arrays
obs_date = obs_time.astype('datetime64[D]')
mod_date = mod_time.astype('datetime64[D]')

# extract indices matching dates by day
print(np.where(np.isin(mod_date, obs_date))[0])
print(np.where(np.isin(mod_date, obs_date)))
matching_indices = np.where(np.isin(mod_date, obs_date))[0]

# Remove mod data do not have matching datetime values to obs
mod_time_2 = mod_time[matching_indices]
mod_temp_2 = mod_temp[matching_indices, :]
mod_salt_2 = mod_salt[matching_indices, :]
mod_date_2 = mod_time_2.astype('datetime64[D]')

#print('no match?')
#matching_indices = np.where(np.isin(mod_date_2, obs_date))[0]
#print(matching_indices.shape)
#inverted = np.logical_not(np.isin(mod_date_2, obs_date))[0]
#print(inverted.shape)
#print('duplicate times?')
#u, c = np.unique(mod_time[matching_indices], return_counts=True)
#dup = u[c > 1]
#print(dup.shape)

print("Shapes below should match now --> ")
print(mod_time_2.shape)
print(obs_time.shape)
print(mod_temp_2.shape)
print(obs_temp.shape)

# There is an undiagnosed problem where it seems there are more modelled data
# than there are observed data - difference is roughly 16 items. 
# can't seem to diagnose this - is it missing data from the server issue? -GO 20240220
# I think this is causing the issues I'm seeing - GO 20240220
# final step to make mod obs match: mod nan where obs nan
# maybe it should go later?
#nan_indices_t = np.isnan(obs_temp)
#nan_indices_s = np.isnan(obs_salt)
#mod_temp_2[nan_indices_t] = np.nan
#mod_salt_2[nan_indices_s] = np.nan

#nan_indices_t = np.isnan(mod_temp_2)
#nan_indices_s = np.isnan(mod_salt_2)
#obs_temp[nan_indices_t] = np.nan
#obs_salt[nan_indices_s] = np.nan


export_plots = False
if export_plots == True:
  print("plotting obs and mod data that should match)")
  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_t = plt.pcolormesh(np.transpose(obs_temp), vmin = 4, vmax = 15, cmap = cm.cm.thermal)
  plt.colorbar(w_t, label = 'deg C')
  axs.invert_yaxis()
  fig.savefig(shared_p + out_p + 'obs_JdF_temp.jpg')

  fact = 0.5
  fig, axs = plt.subplots(1,1, figsize=(27*fact, 8*fact), facecolor='w', edgecolor='k')
  w_t = plt.pcolormesh(np.transpose(mod_temp_2), vmin = 4, vmax = 15, cmap = cm.cm.thermal)
  plt.colorbar(w_t, label = 'deg C')
  axs.invert_yaxis()
  fig.savefig(shared_p + out_p + 'ORAS5_JdF_temp.jpg')


# //////////////////////////////////////////////////////////////////////
# model and obs daily values
# //////////////////////////////////////////////////////////////////////

print("calculating daily mean values for model and obs")

# transform array above so that dates can be plotted with blanks between
# filter by dates corresponding to model minus spin up year
# and bin into daily blocks

start_date1 = '1980-01-01'
end_date1 = '2018-12-31'

dates_mod = pd.DatetimeIndex(mod_time_2)
mod_yrs = dates_mod.year; mod_months = dates_mod.month; mod_days = dates_mod.day

dates_obs = pd.DatetimeIndex(obs_time)
obs_yrs = dates_obs.year; obs_months = dates_obs.month; obs_days = dates_obs.day

dates_all = pd.date_range(start=start_date1, end=end_date1, freq='D')
dates_yrs = dates_all.year; dates_months = dates_all.month; dates_days = dates_all.day

# empty arrays to initialize
mod_salt_ts = np.zeros([len(dates_all), 40]); mod_salt_ts[:] = np.nan
mod_temp_ts = np.zeros([len(dates_all), 40]); mod_temp_ts[:] = np.nan
obs_salt_ts = np.zeros([len(dates_all), 40]); obs_salt_ts[:] = np.nan
obs_temp_ts = np.zeros([len(dates_all), 40]); obs_temp_ts[:] = np.nan

### retreive all model results for a given day. they are given in size [1, x, 40] 
#so you need to take 0 in the first dimension 
for i in range(0, len(dates_all)):

    find_mod = np.where((mod_yrs == dates_yrs[i]) & (mod_months == dates_months[i]) & (mod_days == dates_days[i]))
    find_obs = np.where((obs_yrs == dates_yrs[i]) & (obs_months == dates_months[i]) & (obs_days == dates_days[i]))
    
    
    if (len(dates_mod[find_mod])>0):
         
       mod_available_salt = mod_salt_2[find_mod,:]
       mod_available_salt_toav = mod_available_salt[0,:,:]  
       mod_available_temp = mod_temp_2[find_mod,:]
       mod_available_temp_toav = mod_available_temp[0,:,:]

       obs_available_salt = obs_salt[find_obs,:]
       obs_available_salt_toav = obs_available_salt[0,:,:]
       obs_available_temp = obs_temp[find_obs,:]
       obs_available_temp_toav = obs_available_temp[0,:,:]

       # I expect to see RuntimeWarnings in this block np.nanmean of all nans sometimes
       with warnings.catch_warnings():
         warnings.simplefilter("ignore", category=RuntimeWarning)
         
         mod_salt_ts[i,:] = np.nanmean(mod_available_salt_toav, axis = 0)
         mod_temp_ts[i,:] = np.nanmean(mod_available_temp_toav, axis = 0)
         obs_salt_ts[i,:] = np.nanmean(obs_available_salt_toav, axis = 0)
         obs_temp_ts[i,:] = np.nanmean(obs_available_temp_toav, axis = 0)
           
export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("model and obs temp and salin (4 plots)")
    fig, axs = plt.subplots(4, 1, figsize = (20,16))
    axs = axs.ravel()
    
    w = axs[0].pcolormesh(dates_all, deplevs, np.transpose(mod_temp_ts), vmin=6, vmax=12, cmap = cm.cm.thermal)
    plt.colorbar(w, ax = axs[0])
    axs[0].set_ylim(0, 200)
    axs[0].invert_yaxis()
    
    w = axs[1].pcolormesh(dates_all, deplevs, np.transpose(mod_salt_ts), vmin=28, vmax=31, cmap = cm.cm.haline)
    plt.colorbar(w, ax = axs[1])
    axs[1].set_ylim(0, 200)
    axs[1].invert_yaxis()
    
    w = axs[2].pcolormesh(dates_all, deplevs, np.transpose(obs_temp_ts), vmin=6, vmax=12, cmap = cm.cm.thermal)
    plt.colorbar(w, ax = axs[2])
    axs[2].set_ylim(0, 200)
    axs[2].invert_yaxis()
    
    w = axs[3].pcolormesh(dates_all, deplevs, np.transpose(obs_salt_ts), vmin=28, vmax=31, cmap = cm.cm.haline)
    plt.colorbar(w, ax = axs[3])
    axs[3].set_ylim(0, 200)
    axs[3].invert_yaxis()
    
    axs[0].set_title('Model (ORAS5) Temperature (cons T) at JdF (daily avg)')
    axs[1].set_title('Model (ORAS5) Salinity (abs S, g/kg) at JdF (daily avg)')
    axs[2].set_title('Obs (IOS CTD) Temperature (cons T) at JdF (daily avg)')
    axs[3].set_title('Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (daily avg)')
    
    #plt.tight_layout()
    fig.savefig(shared_p + out_p + 'ORAS5_JdF_temp_salt_withgaps.jpg')

# //////////////////////////////////////////////////////////////////////
# model - obs comparison
# //////////////////////////////////////////////////////////////////////

print("running model - obs comparison")

modobs_temp_ts = np.subtract(mod_temp_ts, obs_temp_ts)
modobs_salt_ts = np.subtract(mod_salt_ts, obs_salt_ts)

export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (daily)")
    fig, axs = plt.subplots(2, 1, figsize = (20,8))
    axs = axs.ravel()
    
    w = axs[0].pcolormesh(dates_all, deplevs, np.transpose(modobs_temp_ts), vmin=-4, vmax=4, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[0])
    axs[0].set_ylim(0, 200)
    axs[0].invert_yaxis()
    
    w = axs[1].pcolormesh(dates_all, deplevs, np.transpose(modobs_salt_ts), vmin=-5, vmax=5, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[1])
    axs[1].set_ylim(0, 200)
    axs[1].invert_yaxis()
    
    axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (daily avg)')
    axs[1].set_title('Model (ORAS5) - Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (daily avg)')
    #plt.tight_layout()
    fig.savefig('Mod-Obs_JdF_temp_salt_withgaps_daily.jpg')

# //////////////////////////////////////////////////////////////////////
# bin using avg monthly (yrs separate), plot
# //////////////////////////////////////////////////////////////////////
print("finding mod-obs avg each month over ts")

dates_all_m = pd.date_range(start=start_date1, end=end_date1, freq='M')
dates_yrs = dates_all_m.year
dates_months = dates_all_m.month

dates_modobs = pd.DatetimeIndex(dates_all)
modobs_yrs = dates_modobs.year 
modobs_months = dates_modobs.month

modobs_temp_ts_mo = np.zeros((dates_all_m.shape[0], modobs_temp_ts.shape[1]))
modobs_salt_ts_mo = np.zeros((dates_all_m.shape[0], modobs_salt_ts.shape[1]))
modobs_temp_ts_mo[:] = np.nan
modobs_salt_ts_mo[:] = np.nan

for i in range(0, len(dates_all_m)):
  find_modobs = np.where((modobs_yrs == dates_yrs[i]) & (modobs_months == dates_months[i]))
  
  if (len(dates_modobs[find_modobs])>0):       
    
    modobs_available_salt = modobs_salt_ts[find_modobs,:]
    modobs_available_salt_toav = modobs_available_salt[0,:,:]  
    
    modobs_available_temp = modobs_temp_ts[find_modobs,:]
    modobs_available_temp_toav = modobs_available_temp[0,:,:]

    # I expect to see RuntimeWarnings in this block np.nanmean of all nans sometimes
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
      
      modobs_salt_ts_mo[i,:] = np.nanmean(modobs_available_salt_toav, axis = 0)
      modobs_temp_ts_mo[i,:] = np.nanmean(modobs_available_temp_toav, axis = 0)

export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (monthly)")
    fig, axs = plt.subplots(2, 1, figsize = (20,8))
    axs = axs.ravel()
    
    w = axs[0].pcolormesh(dates_all_m, deplevs, np.transpose(modobs_temp_ts_mo), vmin=-4, vmax=4, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[0])
    axs[0].set_ylim(0, 200)
    axs[0].invert_yaxis()
    
    w = axs[1].pcolormesh(dates_all_m, deplevs, np.transpose(modobs_salt_ts_mo), vmin=-5, vmax=5, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[1])
    axs[1].set_ylim(0, 200)
    axs[1].invert_yaxis()
    
    axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (monthly avg)')
    axs[1].set_title('Model (ORAS5) - Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (monthly avg)')
    #plt.tight_layout()
    fig.savefig(shared_p + out_p + 'Mod-Obs_JdF_temp_salt_wgaps_monthly.jpg')


# //////////////////////////////////////////////////////////////////////
# bin using avg monthly over all years, plot
# //////////////////////////////////////////////////////////////////////
print("finding mean mod-obs for months across all years.")

mos = np.arange(1,13,1)

dates_modobs = pd.DatetimeIndex(dates_all_m)
modobs_months = dates_modobs.month

modobs_temp_ts_mo_all = np.zeros((mos.shape[0], modobs_temp_ts_mo.shape[1]))
modobs_salt_ts_mo_all = np.zeros((mos.shape[0], modobs_salt_ts_mo.shape[1]))
modobs_temp_ts_mo_all[:] = np.nan
modobs_salt_ts_mo_all[:] = np.nan

for i in range(0, len(mos)):
  find_modobs = np.where(modobs_months == mos[i])
  #print("found this many months with data across all years ", len(dates_modobs[find_modobs]))
  #print("for month ", mos[i])
  
  if (len(dates_modobs[find_modobs])>0):       
    
    modobs_available_salt = modobs_salt_ts_mo[find_modobs,:]
    modobs_available_salt_toav = modobs_available_salt[0,:,:]  
    
    modobs_available_temp = modobs_temp_ts_mo[find_modobs,:]
    modobs_available_temp_toav = modobs_available_temp[0,:,:]

    # I expect to see RuntimeWarnings in this block np.nanmean of all nans sometimes
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
      
      modobs_salt_ts_mo_all[i,:] = np.nanmean(modobs_available_salt_toav, axis = 0)
      modobs_temp_ts_mo_all[i,:] = np.nanmean(modobs_available_temp_toav, axis = 0)


export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (monthly over all yrs)")
    fig, axs = plt.subplots(2, 1, figsize = (6,6))
    axs = axs.ravel()
    
    my_xticks_lab = ['J','F','M','A','M','J','J','A','S','O','N','D']
    
    w = axs[0].pcolormesh(mos, deplevs, np.transpose(modobs_temp_ts_mo_all), vmin=-4, vmax=4, cmap = cm.cm.balance)
    cb1 = plt.colorbar(w, ax = axs[0])
    cb1.ax.set_ylabel('Temp ($^\circ$C)', fontdict={'fontsize':10, 'rotation':270}, labelpad=10)
    axs[0].set_ylim(0, 175)
    axs[0].invert_yaxis()
    axs[0].set_ylabel('Depth (m)', fontsize=10)
    axs[0].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    axs[0].set_xticklabels(my_xticks_lab)
    
    w = axs[1].pcolormesh(mos, deplevs, np.transpose(modobs_salt_ts_mo_all), vmin=-5, vmax=5, cmap = cm.cm.balance)
    cb2 = plt.colorbar(w, ax = axs[1])
    cb2.ax.set_ylabel('Salinity (PSU)', fontdict={'fontsize':10, 'rotation':270}, labelpad=10)
    axs[1].set_ylim(0, 175)
    axs[1].invert_yaxis()
    axs[1].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    axs[1].set_xticklabels(my_xticks_lab)
    print(axs[1].get_xticks())
    axs[1].set_xlabel("Month", fontsize=10)
    axs[1].set_ylabel('Depth (m)', fontsize=10)
    
    #axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (monthly avg across yrs)')
    axs[0].set_title('Temperature Bias (ORAS5)', fontsize=11)
    #axs[1].set_title('Model (ORAS5) - Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (monthly avg across yrs)')
    axs[1].set_title('Salinity Bias (ORAS5)', fontsize=11)
    plt.tight_layout()
    fig.savefig(shared_p + out_p + 'Mod-Obs_JdF_temp_salt_monthly_allyrs.jpg')
    
export_plots = True
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (monthly over all yrs)")
    fig, ax = plt.subplots(1, 1, figsize = (6,3)) #w,h
    
    my_xticks_lab = ['J','F','M','A','M','J','J','A','S','O','N','D']
    
    w = ax.pcolormesh(mos, deplevs, np.transpose(modobs_temp_ts_mo_all), vmin=-4, vmax=4, cmap = cm.cm.balance)
    cb1 = plt.colorbar(w, ax = ax)
    cb1.ax.set_ylabel('Temp ($^\circ$C)', fontdict={'fontsize':10, 'rotation':270}, labelpad=10)
    ax.set_ylim(0, 175)
    ax.invert_yaxis()
    ax.set_ylabel('Depth (m)', fontsize=10)
    ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ax.set_xticklabels(my_xticks_lab)
    
    #axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (monthly avg across yrs)')
    #ax.set_title('Temperature Bias (ORAS5)', fontsize=11)
    ax.set_xlabel("Month", fontsize=10)
    plt.tight_layout()
    fig.savefig(shared_p + out_p + 'Fig03-Mod-Obs_JdF_temp_monthly_allyrs.png')
    fig.savefig(shared_p + out_p + 'Fig03-Mod-Obs_JdF_temp_monthly_allyrs.pdf')
    fig.savefig(shared_p + out_p + 'Fig03-Mod-Obs_JdF_temp_monthly_allyrs.svg')

    

    

# to-do: export .nc file    
    
# save as netcdf
savenam = shared_p + out_p + 'JdF_mod-obs_ORAS_monthlyallyrs.nc'
data_vars = {
    'salinity_bias':(['time_counter', 'deptht'], modobs_salt_ts_mo_all,
    {'units': 'psu',
    'long_name':'practical salinity'}),
    'temperature_bias':(['time_counter', 'deptht'], modobs_temp_ts_mo_all,
    {'units': 'deg C',
    'long_name':'in situ temperature'}),
    'month':(['time_counter'], mos,
    {'units': '1-12',
    'long_name':''}),             
    }

# define coordinates
coords = {'time_counter': (['time_counter'], mos),
    'deptht': (['deptht'], deplevs)}
# define global attributes
attrs = {'made in':'ORAS_4_processObs.py',
'desc': 'monthly avg ORAS model minus CTD obs (IOS) at mouth of JdF, 1980 - 2018'
}
ds = xr.Dataset(data_vars=data_vars, coords=coords,attrs=attrs)
ds.to_netcdf(savenam)

# //////////////////////////////////////////////////////////////////////
# bin using avg monthly over all two time periods (1980 - 1999, 2000 - 2018), plot
# //////////////////////////////////////////////////////////////////////

modobs_yrs = dates_modobs.year

#Same as above grouped pre 2000 and post 2000
modobs_temp_ts_mo_all = np.zeros((mos.shape[0], modobs_temp_ts_mo.shape[1]))
modobs_salt_ts_mo_all = np.zeros((mos.shape[0], modobs_salt_ts_mo.shape[1]))
modobs_temp_ts_mo_all[:] = np.nan
modobs_salt_ts_mo_all[:] = np.nan

for i in range(0, len(mos)):
  find_modobs = np.where((modobs_months == mos[i]) & (modobs_yrs < 2000))
  #print("found this many months with data across all years ", len(dates_modobs[find_modobs]))
  #print("for month ", mos[i])
  
  if (len(dates_modobs[find_modobs])>0):       
    
    modobs_available_salt = modobs_salt_ts_mo[find_modobs,:]
    modobs_available_salt_toav = modobs_available_salt[0,:,:]  
    
    modobs_available_temp = modobs_temp_ts_mo[find_modobs,:]
    modobs_available_temp_toav = modobs_available_temp[0,:,:]

    # I expect to see RuntimeWarnings in this block np.nanmean of all nans sometimes
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
      
      modobs_salt_ts_mo_all[i,:] = np.nanmean(modobs_available_salt_toav, axis = 0)
      modobs_temp_ts_mo_all[i,:] = np.nanmean(modobs_available_temp_toav, axis = 0)


export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (monthly over 1980-1999)")
    fig, axs = plt.subplots(2, 1, figsize = (10,8))
    axs = axs.ravel()
    
    w = axs[0].pcolormesh(mos, deplevs, np.transpose(modobs_temp_ts_mo_all), vmin=-4, vmax=4, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[0])
    axs[0].set_ylim(0, 200)
    axs[0].invert_yaxis()
    
    w = axs[1].pcolormesh(mos, deplevs, np.transpose(modobs_salt_ts_mo_all), vmin=-5, vmax=5, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[1])
    axs[1].set_ylim(0, 200)
    axs[1].invert_yaxis()
    
    axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (monthly avg across 1980 - 1999)')
    axs[1].set_title('Model (ORAS5) - Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (monthly avg across 1980 - 1999)')
    #plt.tight_layout()
    fig.savefig(shared_p + out_p + 'Mod-Obs_JdF_temp_salt_monthly_1980to1999.jpg')
    
modobs_temp_ts_mo_all = np.zeros((mos.shape[0], modobs_temp_ts_mo.shape[1]))
modobs_salt_ts_mo_all = np.zeros((mos.shape[0], modobs_salt_ts_mo.shape[1]))
modobs_temp_ts_mo_all[:] = np.nan
modobs_salt_ts_mo_all[:] = np.nan

for i in range(0, len(mos)):
  find_modobs = np.where((modobs_months == mos[i]) & (modobs_yrs >= 2000))
  #print("found this many months with data across all years ", len(dates_modobs[find_modobs]))
  #print("for month ", mos[i])
  
  if (len(dates_modobs[find_modobs])>0):       
    
    modobs_available_salt = modobs_salt_ts_mo[find_modobs,:]
    modobs_available_salt_toav = modobs_available_salt[0,:,:]  
    
    modobs_available_temp = modobs_temp_ts_mo[find_modobs,:]
    modobs_available_temp_toav = modobs_available_temp[0,:,:]

    # I expect to see RuntimeWarnings in this block np.nanmean of all nans sometimes
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
      
      modobs_salt_ts_mo_all[i,:] = np.nanmean(modobs_available_salt_toav, axis = 0)
      modobs_temp_ts_mo_all[i,:] = np.nanmean(modobs_available_temp_toav, axis = 0)


export_plots = False
if export_plots == True:
    # visualize, include gaps
    print("plotting model - obs (monthly over 2000-2018)")
    fig, axs = plt.subplots(2, 1, figsize = (10,8))
    axs = axs.ravel()
    
    w = axs[0].pcolormesh(mos, deplevs, np.transpose(modobs_temp_ts_mo_all), vmin=-4, vmax=4, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[0])
    axs[0].set_ylim(0, 200)
    axs[0].invert_yaxis()
    
    w = axs[1].pcolormesh(mos, deplevs, np.transpose(modobs_salt_ts_mo_all), vmin=-5, vmax=5, cmap = cm.cm.balance)
    plt.colorbar(w, ax = axs[1])
    axs[1].set_ylim(0, 200)
    axs[1].invert_yaxis()
    
    axs[0].set_title('Model (ORAS5) - Obs (IOS CTD) Temperature (cons T) at JdF (monthly avg across 2000 - 2018)')
    axs[1].set_title('Model (ORAS5) - Obs (IOS CTD) Salinity (abs S, g/kg) at JdF (monthly avg across 2000 - 2018)')
    #plt.tight_layout()
    fig.savefig(shared_p + out_p + 'Mod-Obs_JdF_temp_salt_monthly_2000to2018.jpg') 
    
print("Finished.") 


#