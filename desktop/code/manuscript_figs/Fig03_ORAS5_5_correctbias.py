# adjusts boundary files that are already prepped for nemo
# using a bias correction factor, precalculated (monthly over depths).
# applies same correction across all boundary cells
# - G Oldford 2023-04-27
# - last edited: 2023-07-12
# this script meant to be run on unix (server: graham)
# downloaded 2024-03-05


import netCDF4 as nc
import numpy as np
import os
import shutil
import sys

# files look like 'ts_W_y{}m{}.nc'.format(yr,mo) and are shape(1,40,10,32) (time,z,y,x)
bdy_in_p = '/project/6006412/mdunphy/Forcing/forcing_oras5_20210501/bdy_ts/'
bdy_out_p = '/project/6006412/goldford/data/Forcing/oras5_bdy_ts_biasfix_t_jul17/'

shared_p = '/project/6006412/goldford/'
clim_p = 'ANALYSIS/ORAS5-JdF/CUSTOM-GO/'
clim_f = 'JdF_mod-obs_ORAS_monthlyallyrs.nc'

bias_correct_temp = True
bias_correct_salt = False

var_clim_t = 'temperature_bias'
var_clim_s = 'salinity_bias'
var_bdy_t = 'votemper'
var_bdy_s = 'vosaline'

# get bias climatology from nc
fix_dep1 = True # adjustment opt, no bias correction where z=1
fix_dep29up = True # adjustment opt, add bias correction around 147 m (error in creating clim file) -20230712
with nc.Dataset(os.path.join(shared_p, clim_p, clim_f)) as clim_d:

  clim_d_t = np.zeros(clim_d[var_clim_t].shape)
  clim_d_t[:,:] = np.nan
  clim_d_t[:,:] = clim_d[var_clim_t][:,:]
  
  clim_d_s = np.zeros(clim_d[var_clim_s].shape)
  clim_d_s[:,:] = np.nan
  clim_d_s[:,:] = clim_d[var_clim_s][:,:]
  
  # optional adjustment
  if fix_dep1 == True:
    # if first depth always np.nan, use bias correction from z=2
    clim_d_t[:,0] = clim_d_t[:,1]
    clim_d_s[:,0] = clim_d_s[:,1]
    
  if fix_dep29up == True:
    # no clima bias adjust for z 29, copy z 28 (temp hack fix 20230712)
    #clim_d_t[:,29] = clim_d_t[:,28]
    clim_d_t[:,29:] = clim_d_t[:, 28][:, None] 
    clim_d_s[:,29:] = clim_d_s[:, 28][:, None] 


for yr in np.arange(1979, 2020, 1):
  print("processing year ", yr)
  mo_n = 0
  for mo in ["01","02","03","04","05","06","07","08","09","10","11","12"]:
  
    # input files have only one 2019 file for 2019-01 - this catches to stop before 2019-02
    if yr == 2019 and mo_n >= 1:
      break
  
    # get data for month from climatology mod-obs
    clim_mo_t = clim_d_t[mo_n,:]
    clim_mo_s = clim_d_s[mo_n,:]
    
    clim_mo_t_reshp = clim_mo_t.reshape(1, 40, 1, 1)
    clim_mo_s_reshp = clim_mo_s.reshape(1, 40, 1, 1)
    
    infile = os.path.join(bdy_in_p, 'ts_W_y{}m{}.nc'.format(yr,mo))
    
    if not os.path.exists(infile):  print("No input file (exiting). Not found: " + infile); sys.exit()
    
    outfile = os.path.join(bdy_out_p, 'ts_W_y{}m{}_adjst_t.nc'.format(yr,mo))
    shutil.copy2(infile,outfile)  # copy original file

    with nc.Dataset(outfile, 'r+') as ncf:
      #1x40x10x32 arrays (time, z, y, x) for t,s
      
      # temp
      #bdy_t = np.zeros(ncf[var_bdy_t].shape)
      #bdy_t[:,:,:,:] = ncf[var_bdy_t][:,:,:,:]
      
      bdy_t = ncf[var_bdy_t][:,:,:,:]
      
      if bias_correct_temp == True:
        bdy_t[bdy_t==0] = np.nan # temporarily switch to nan
        bdy_t_new = bdy_t - clim_mo_t_reshp
        bdy_t_new[np.isnan(bdy_t_new)] = 0
      else: 
        bdy_t_new = bdy_t
      
      # salt
      #bdy_s = np.zeros(ncf[var_bdy_s].shape)
      #bdy_s[:,:,:,:] = ncf[var_bdy_s][:,:,:,:]
      bdy_s = ncf[var_bdy_s][:,:,:,:]
      
      if bias_correct_salt == True:
        bdy_s[bdy_s==0] = np.nan
        bdy_s_new = bdy_s - clim_mo_s_reshp
        bdy_s_new[np.isnan(bdy_s_new)] = 0
      else: 
        bdy_s_new = bdy_s
      
      # write bias-corrected data
      ncf[var_bdy_t][:,:,:,:] = bdy_t_new
      ncf[var_bdy_s][:,:,:,:] = bdy_s_new
      
    mo_n += 1

# debug    
#    if mo_n == 1:
#      break



#print("temp BEFORE:")
#print(bdy_t[0,:,4,1])
#print("temp bias:")
#print(clim_mo_t_reshp)
#print("temp AFTER:")
#print(bdy_t_new[0,:,4,1])
#
#print("salt BEFORE:")
#print(bdy_s[0,:,4,1])
#print("salt bias:")
#print(clim_mo_s_reshp)
#print("salt AFTER:")
#print(bdy_s_new[0,:,4,1])
