import numpy as np
import netCDF4 as nc
import os
from matplotlib import pyplot as plt

path_wi='/project/6006412/goldford/ECOSPACE/DATA/wind_hourly/'
wind_f='RDRS21_NEMOgrid_wind_1999.nc'

path_ma='/project/6006412/goldford/data/grid/'
mask_f='mesh_mask_20210406.nc'

with nc.Dataset(os.path.join(path_ma,mask_f)) as ma_nc:
  #nlat = ma_nc.variables['gphit'][0,...]   # neglected the 1/2 grid box shift here
  #nlon = ma_nc.variables['glamt'][0,...]
  tmask = ma_nc.variables["tmask"][0,0,:, :].filled()
 
  
print(tmask.shape)


with nc.Dataset(os.path.join(path_wi,wind_f)) as wi_nc:
  #nlat = li_nc.variables['gphit'][0,...]   # neglected the 1/2 grid box shift here
  #nlon = li_nc.variables['glamt'][0,...]
  v10m = wi_nc.variables["v10m"][1000,...]
  u10m = wi_nc.variables["u10m"][1000,...]
  
print(v10m.shape)


fig, ax = plt.subplots(figsize=(8, 5))
#masked_wi = np.multiply(v10m[:,:],tmask[:,:])

#w = ax.pcolormesh(v10m[:,:], cmap=plt.pink())
#cb = fig.colorbar(w, ax = ax,label='[units]')

print(type(v10m))

v10m_neg = np.where(v10m > 0, 1, 0)
#v10m_neg = np.where(v10m_neg > 3, 2, v10m_neg)
#v10m_neg = np.where(v10m_neg > 6, 3, v10m_neg)

#v10m_neg_ma = np.ma.masked_array(v10m_neg, v10m_neg == 1)
w3 = ax.pcolormesh(v10m_neg[:,:], cmap=plt.pink())
cb3 = fig.colorbar(w3, ax = ax,label='[units]')

#tmask_invert = (tmask - 1) * -1
#tmask_inv_ma = np.ma.masked_array(tmask, tmask == 1)
#w2 = ax.pcolormesh(tmask_inv_ma[:,:], cmap=plt.cm.binary)
#cb2 = fig.colorbar(w2, ax = ax,label='mask')


plt.tight_layout()
plt.title("lwind_check")
plt.savefig("windcheck3.png")
