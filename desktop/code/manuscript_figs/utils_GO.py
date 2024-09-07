# created 2023-08-17 by GO
# keep helper code functions etc out of notebooks
from datetime import datetime, timedelta
import numpy as np
import calendar
from stats_GO import willmott1981
import skill_metrics as sm
import os
import netCDF4 as nc
from scipy import interpolate
import matplotlib.pyplot as plt
import xarray as xr


class new_utils_GO():

  # interpolated to given depths, takes nc as input
  # adapted from fn in Nanoose ipynb's 202304  
  def get_model_interpolated(ncname, mod_dlevs, verbose = False):
  
      tobs = xr.open_dataset(ncname)
      obs_d = tobs['Pres'].values
      obs_t = tobs['cTemp'].values
      obs_s = tobs['aSal'].values
      ttime = tobs['time'][0].values
      
      #throw out stuff below 400
      filt = obs_d > 400
      obs_d[filt] = np.nan
      obs_t[filt] = np.nan
      obs_s[filt] = np.nan
      
      try:
          f = interpolate.interp1d(obs_d, obs_t) #temperature
          f2 = interpolate.interp1d(obs_d, obs_s) #salinity
  
          ## can only interpolate to model points that are within observations 
          mod_d = mod_dlevs[(mod_dlevs<max(obs_d)) & (mod_dlevs>min(obs_d))]
          firstind = np.where(mod_dlevs == np.min(mod_d))[0][0] ## first model index we were able to interpolate to 
          interp_t = f(mod_d)   # use interpolation function returned by `interp1d`
          interp_s = f2(mod_d)
  
          ### 
          t_full = np.zeros(40)
          t_full[:] = -999
          t_full[firstind:firstind+len(mod_d)] = interp_t
          t_full[t_full < -900] = np.nan
          t_full[t_full > 100] = np.nan
  
          s_full = np.zeros(40)
          s_full[:] = -999
          s_full[firstind:firstind+len(mod_d)] = interp_s
          s_full[s_full < -900] = np.nan
          s_full[s_full > 100] = np.nan
  
          if verbose:
              fig, axs =  plt.subplots(1,2)
              axs = axs.ravel()
              axs[0].plot(obs_t, obs_d, 'ob', interp_t, mod_d, 'or')
              axs[1].plot(obs_s, obs_d, 'ob', interp_s, mod_d, 'or')
              axs[0].set_title('temperature cons')
              axs[1].set_title('salinity abs')
              axs[0].invert_yaxis()
              axs[1].invert_yaxis()
  
              fig.suptitle(f'blue is CTD, red is interp. to mod depths \n date of cast: {ttime}')
              plt.tight_layout()
              plt.show()
      except:
                 
          t_full = np.zeros(40)
          t_full[:] = -999
          t_full[t_full < -900] = np.nan
          s_full = np.zeros(40)
          s_full[:] = -999
          s_full[s_full < -900] = np.nan
          
      return t_full, s_full, ttime # not sure ttime is needed -GO 202304
  
  # interpolated to given depths, takes np arrays as input
  # adapted from fn in Nanoose ipynb's 202304  
  def get_model_interpolated_ar(d_pres, d_salt, d_temp, d_time, mod_dlevs, verbose = False):
  
      
      obs_d = d_pres
      obs_t = d_temp
      obs_s = d_salt
      ttime = d_time
      
      #throw out stuff below 400
      filt = obs_d > 400
      obs_d[filt] = np.nan
      obs_t[filt] = np.nan
      obs_s[filt] = np.nan
      
      
      try:
          f = interpolate.interp1d(obs_d, obs_t) #temperature
          f2 = interpolate.interp1d(obs_d, obs_s) #salinity
  
          ## can only interpolate to model points that are within observations 
          mod_d = mod_dlevs[(mod_dlevs<max(obs_d)) & (mod_dlevs>min(obs_d))]
          firstind = np.where(mod_dlevs == np.min(mod_d))[0][0] ## first model index we were able to interpolate to 
          interp_t = f(mod_d)   # use interpolation function returned by `interp1d`
          interp_s = f2(mod_d)
  
          ### 
          t_full = np.zeros(40)
          t_full[:] = -999
          t_full[firstind:firstind+len(mod_d)] = interp_t
          t_full[t_full < -900] = np.nan
          t_full[t_full > 100] = np.nan
  
          s_full = np.zeros(40)
          s_full[:] = -999
          s_full[firstind:firstind+len(mod_d)] = interp_s
          s_full[s_full < -900] = np.nan
          s_full[s_full > 100] = np.nan
  
          if verbose:
              fig, axs =  plt.subplots(1,2)
              axs = axs.ravel()
              axs[0].plot(obs_t, obs_d, 'ob', interp_t, mod_d, 'or')
              axs[1].plot(obs_s, obs_d, 'ob', interp_s, mod_d, 'or')
              axs[0].set_title('temperature cons')
              axs[1].set_title('salinity abs')
              axs[0].invert_yaxis()
              axs[1].invert_yaxis()
  
              fig.suptitle(f'blue is CTD, red is interp. to mod depths \n date of cast: {ttime}')
              plt.tight_layout()
              plt.show()
      except:
                 
          t_full = np.zeros(40)
          t_full[:] = -999
          t_full[t_full < -900] = np.nan
          s_full = np.zeros(40)
          s_full[:] = -999
          s_full[s_full < -900] = np.nan
          
      return t_full, s_full
          
# class write_to_netcdf from data_conversion_functions - GO 2022-12
class write_to_netcdf_GO(): 

    def CTD_Cast(self, saveDir, station, startT, P, T, consTemp, S, absSal, pDen, dep, lat, lon):
        writeName = station.replace('/', '-') + '_CastCTD_' + str(startT)[:10] + '_' + str(startT)[11:13] + str(startT)[
                                                                                                            14:16] + 'h.nc'

        svPath = general_utils().makePath(os.path.join(saveDir, writeName))

        with nc.Dataset(svPath, 'w', format='NETCDF4') as f:
            f.createDimension('pressure', len(P))
            f.createDimension('single', 1)

            pres = f.createVariable('Pres', 'f4', 'pressure', fill_value=-999)
            pres[:] = P
            pres.units = 'dBar'
            pres.description = 'recorded pressure levels'

            Temp = f.createVariable('temperature', 'f8', 'pressure', fill_value=-999)
            Temp[:] = T
            Temp.type = 'in-situ'
            Temp.units = 'degC'
            Temp.description = 'recorded temperature'

            cT = f.createVariable('cTemp', 'f8', 'pressure', fill_value=-999)
            cT[:] = consTemp
            cT.units = 'degC'
            cT.description = 'conservative temperature'

            Sal = f.createVariable('salinity', 'f8', 'pressure', fill_value=-999)
            Sal[:] = S
            Sal.type = 'practical'
            Sal.units = 'psu'
            Sal.description = 'recorded salinity'

            aSal = f.createVariable('aSal', 'f8', 'pressure', fill_value=-999)
            aSal[:] = absSal
            aSal.units = 'g/kg'
            aSal.description = 'absolute salinity'

            pden = f.createVariable('pDen', 'f8', 'pressure', fill_value=-999)
            pden[:] = pDen
            pden.units = 'kg/m3'
            pden.description = 'potential density'

            time = f.createVariable('time', 'i8', 'single', fill_value=-999)
            time[:] = 0.
            time.units = 'seconds since ' + str(startT)

            d = f.createVariable('depth', 'f4', 'single', fill_value=-999)
            if dep == 0.:
                d[:] = -999.
            else:
                d[:] = dep
            d.units = 'm'
            d.description = 'water depth at cast location'

            Lat = f.createVariable('latitude', 'f4', 'single')
            Lat[:] = lat
            Lat.units = 'degrees north'
            Lat.description = 'cast latitude'

            Lon = f.createVariable('longitude', 'f4', 'single')
            Lon[:] = lon
            Lon.units = 'degrees east'
            Lon.description = 'cast longitude (-180 -> 180 degrees)'

        return svPath
    
    def nanoose_CTD_Cast(self, saveDir, station, startT, P, T, consTemp, S, absSal, pDen, dep, lat, lon):
        
        writeName = station.replace('/', '-') + '_CastCTD_' + str(startT)[:10] + '_' + str(startT)[11:13] + str(startT)[14:16] + 'h.nc'
        svPath = general_utils().makePath(os.path.join(saveDir, writeName))
        
        # GO - 2022-12
        newT = (startT-datetime.datetime(1900,1,1)).total_seconds()
        
        with nc.Dataset(svPath, 'w', format='NETCDF4') as f:
            
            f.createDimension('time_dim', None)
            f.createDimension('lon_dim',1)
            f.createDimension('lat_dim',1)
            f.createDimension('pres_dim', len(P))
            #f.createDimension('single', 1)
    
            # variables
            time_var = f.createVariable('time_var', 'float32', ('time_dim',))
            time_var.units = 'seconds since 1900-01-01 00:00:00'
            #time_var[:] = np.full((len(P)), newT)
            time_var[:] = newT
            
            Lon = f.createVariable('longitude', 'f4', ('lat_dim',))
            Lon[:] = lon
            Lon.units = 'degrees east'
            Lon.description = 'cast longitude (-180 -> 180 degrees)'
            
            Lat = f.createVariable('latitude', 'f4', ('lat_dim',))
            Lat[:] = lat
            Lat.units = 'degrees north'
            Lat.description = 'cast latitude'
    
            # time_counter
#            time_counter = f.createVariable('time', 'float32', ('time_counter'))
#            time_counter.units = 'seconds since 1900-01-01 00:00:00'
#            time_counter[:] = nc.date2num(dates_, time_counter.units)
    
            pres = f.createVariable('pres_var', 'f4', ('pres_dim',))
            pres[:] = P
            pres.units = 'dBar'
            pres.description = 'recorded pressure levels'
    
            Temp = f.createVariable('temperature', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Temp[0,:,0,0] = T
            Temp.type = 'in-situ'
            Temp.units = 'degC'
            Temp.description = 'recorded temperature'
    
            cT = f.createVariable('cTemp', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            cT[0,:,0,0] = consTemp
            cT.units = 'degC'
            cT.description = 'conservative temperature'
    
            Sal = f.createVariable('salinity', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Sal[0,:,0,0] = S
            Sal.type = 'practical'
            Sal.units = 'psu'
            Sal.description = 'recorded salinity'
    
            aSal = f.createVariable('aSal', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            aSal[0,:,0,0] = absSal
            aSal.units = 'g/kg'
            aSal.description = 'absolute salinity'
    
            pden = f.createVariable('pDen', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            pden[0,:,0,0] = pDen
            pden.units = 'kg/m3'
            pden.description = 'potential density'
    
            d = f.createVariable('depth', 'f4', ('time_dim','pres_dim','lon_dim','lat_dim',))
            if dep == 0.:
                d[0,:,0,0] = -999.
            else:
                d[0,:,0,0] = dep
            d.units = 'm'
            d.description = 'water depth at cast location'

        return svPath

    
    # write interpolated nanoose station data (not all are simple CTDs)
    def write_interp_nanoose(self, saveDir, outfile, startT, P, T, consTemp, S, absSal, pDen, dep, lat, lon):
        
        #writeName = station.replace('/', '-') + '_CastCTD_' + str(startT)[:10] + '_' + str(startT)[11:13] + str(startT)[14:16] + 'h.nc'
        writeName = outfile

        #svPath = utils().makePath(os.path.join(saveDir, writeName))
        os.makedirs(saveDir, exist_ok=True)
        fileout = os.path.join(saveDir, outfile)
        
        # GO - 2022-12
        #newT = (startT-datetime.datetime(1900,1,1)).total_seconds()
        newT = startT
        
        # TEMP - testing if ncrcat problem is issue with pressure dimension
        with nc.Dataset(fileout, 'w', format='NETCDF4') as f:
            
            f.createDimension('time_dim', None)
            f.createDimension('lon_dim',1)
            f.createDimension('lat_dim',1)
            f.createDimension('pres_dim', len(P))
            #f.createDimension('single', 1)
    
            # variables
            time_var = f.createVariable('time_var', 'float32', ('time_dim',))
            time_var.units = 'seconds since 1900-01-01 00:00:00'
            #time_var[:] = np.full((len(P)), newT)
            time_var[:] = newT
            
            Lon = f.createVariable('longitude', 'f4', ('lat_dim',))
            Lon[:] = lon
            Lon.units = 'degrees east'
            Lon.description = 'cast longitude (-180 -> 180 degrees)'
            
            Lat = f.createVariable('latitude', 'f4', ('lat_dim',))
            Lat[:] = lat
            Lat.units = 'degrees north'
            Lat.description = 'cast latitude'
    
            # time_counter
#            time_counter = f.createVariable('time', 'float32', ('time_counter'))
#            time_counter.units = 'seconds since 1900-01-01 00:00:00'
#            time_counter[:] = nc.date2num(dates_, time_counter.units)
    
            pres = f.createVariable('pres_var', 'f4', ('pres_dim',))
            pres[:] = P
            pres.units = 'dBar'
            pres.description = 'recorded pressure levels'
    
            Temp = f.createVariable('temperature', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Temp[0,:,0,0] = T
            Temp.type = 'in-situ'
            Temp.units = 'degC'
            Temp.description = 'recorded temperature'
    
            cT = f.createVariable('cTemp', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            cT[0,:,0,0] = consTemp
            cT.units = 'degC'
            cT.description = 'conservative temperature'
    
            Sal = f.createVariable('salinity', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Sal[0,:,0,0] = S
            Sal.type = 'practical'
            Sal.units = 'psu'
            Sal.description = 'recorded salinity'
    
            aSal = f.createVariable('aSal', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            aSal[0,:,0,0] = absSal
            aSal.units = 'g/kg'
            aSal.description = 'absolute salinity'
    
            pden = f.createVariable('pDen', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            pden[0,:,0,0] = pDen
            pden.units = 'kg/m3'
            pden.description = 'potential density'
    
            d = f.createVariable('depth', 'f4', ('time_dim','pres_dim','lon_dim','lat_dim',))
            if dep.all() == 0.:
                d[0,:,0,0] = -999.
            else:
                d[0,:,0,0] = dep
            d.units = 'm'
            d.description = 'water depth at cast location'

        return fileout
  
    # write interpolated nanoose station data (not all are simple CTDs)
    def write_empty_annual(self, saveDir, outfile, startT, P, dep):
        
        writeName = outfile

        os.makedirs(saveDir, exist_ok=True)
        fileout = os.path.join(saveDir, outfile)
        
        # GO - 2022-12
        #newT = (startT-datetime.datetime(1900,1,1)).total_seconds()
        newT = startT
        
        with nc.Dataset(fileout, 'w', format='NETCDF4') as f:
            
            f.createDimension('time_dim', None)
            f.createDimension('lon_dim',1)
            f.createDimension('lat_dim',1)
            f.createDimension('pres_dim', len(P))
            #f.createDimension('single', 1)
    
            # variables
            time_var = f.createVariable('time_var', 'float32', ('time_dim',))
            time_var.units = 'seconds since 1900-01-01 00:00:00'
            #time_var[:] = np.full((len(P)), newT)
            time_var[:] = newT
            
            Lon = f.createVariable('longitude', 'f4', ('lat_dim',), fill_value=-999)
            Lon[:] = -999.
            Lon.units = 'degrees east'
            Lon.description = 'cast longitude (-180 -> 180 degrees)'
            
            Lat = f.createVariable('latitude', 'f4', ('lat_dim',), fill_value=-999)
            Lat[:] = -999.
            Lat.units = 'degrees north'
            Lat.description = 'cast latitude'
    
            # time_counter
#            time_counter = f.createVariable('time', 'float32', ('time_counter'))
#            time_counter.units = 'seconds since 1900-01-01 00:00:00'
#            time_counter[:] = nc.date2num(dates_, time_counter.units)
    
            pres = f.createVariable('pres_var', 'f4', ('pres_dim',), fill_value=-999)
            pres[:] = P
            pres.units = 'dBar'
            pres.description = 'recorded pressure levels'
    
            Temp = f.createVariable('temperature', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Temp[:] = -999.
            Temp.type = 'in-situ'
            Temp.units = 'degC'
            Temp.description = 'recorded temperature'
    
            cT = f.createVariable('cTemp', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            cT[:] = -999.
            cT.units = 'degC'
            cT.description = 'conservative temperature'
    
            Sal = f.createVariable('salinity', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            Sal[:] = -999.
            Sal.type = 'practical'
            Sal.units = 'psu'
            Sal.description = 'recorded salinity'
    
            aSal = f.createVariable('aSal', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            aSal[:] = -999.
            aSal.units = 'g/kg'
            aSal.description = 'absolute salinity'
    
            pden = f.createVariable('pDen', 'f8', ('time_dim','pres_dim','lon_dim','lat_dim',), fill_value=-999)
            pden[:] = -999.
            pden.units = 'kg/m3'
            pden.description = 'potential density'
    
            d = f.createVariable('depth', 'f4', ('pres_dim',), fill_value=-999)
            if dep.all() == 0.:
                d[:] = -999.
            else:
                d[:] = dep
            d.units = 'm'
            d.description = 'water depth at cast location'

        return fileout
 
  
  
  #Based on code in https://gitlab.com/FA12/python-analysis-package/-/blob/master/analysispkg/pkg_interp.py
  #purpose: vertically interpolate profiles to nemo model vert levels
  # from pkg_interp.py
class pkg_interp_GO(): 

       # example code from pyap (not using this)
       # https://gitlab.com/FA12/python-analysis-package/-/blob/master/analysispkg/pkg_extract_forecast.py
       #u = u.reshape(nz, -1)
       #v = v.reshape(nz, -1)
       #i1, i2, w1, w2 = pkg_interp.interp_weights_1d(zi_tiled, zut, zmask=None, extrap_threshold=0)
       #u = np.take(u, i1) * w1 + np.take(u, i2) * w2  # where u has z.shape
       #u = u.reshape(nx, ntime)
       
       # open model file to get levels
       # for each obs file, open, interp, re-wrte
       # levels

  def interp_weights_1d(zi, z, zmask=None, extrap_threshold=None):
      """ 1d linear interpolation indices and weights.
  
      Works on n-d arrays along 1st dimension.
  
      Parameters
      ----------
          zi : array_like of shape (n,d1,...)
              N-d array or scalar. Coordinate of `n` points to interpolate to.
          z : array_like of shape (m,d1,...)
              Points to interpolate from. Second and subsequent dimensions must
              match those of `zi`.
          zmask : array_like of shape (m,d1,...), optional
              Mask for `z` with 0's for invalid points.
          extrap_threshold : float, optional
              Extrapolate no further that the threshold (same units as `z` and `zi`).
              No threshold by default.
  
      Returns
      -------
          i1 : array_like of shape (n,d1,...)
          i2 : array_like of shape (n,d1,...)
              Interpolation indices, ravelled to use with (m,d1,...) arrays.
          w1 : array_like of shape (n,d1,...)
          w2 : array_like of shape (n,d1,...)
              Corresponding weights.
  
      Example
      -------
      To apply indices and weights:
          vi = np.take(v,i1)*w1 + np.take(v,i2)*w2 # where v has z.shape
      """
      # generalize for inputs of various shapes, including scalars
      scalar = np.isscalar(zi)
      if not scalar:
          sz = zi.shape
      zi, z = atleast_2d0(zi, z)
  
      n = zi.shape[0]
      i1, i2 = [np.zeros(zi.shape, dtype=np.int32) for i in range(2)]  # initialize
      w1, w2 = [np.full(zi.shape, np.nan) for i in range(2)]  # initialize
  
      # deal with completely masked nodes
      if zmask is not None:
          allmasked = np.all(zmask == 0, axis=0)
          nodes = np.where(~allmasked)[0]  # not masked
      else:
          nodes = range(zi.shape[1])  # all nodes
  
      for kl in nodes:
          if zmask is not None:
              zm = zmask[:, kl]
              i1[:, kl], i2[:, kl], w1[:, kl], w2[:, kl] = vinterp1d(zi[:, kl], z[zm, kl], extrap_threshold)
          else:
              i1[:, kl], i2[:, kl], w1[:, kl], w2[:, kl] = vinterp1d(zi[:, kl], z[:, kl], extrap_threshold)
  
      # ravel indices for subsequent indexing of arrays
      dim1 = zi.shape[1:]  # 2nd and subsequent dimensions
      dim1prod = np.prod(dim1)  # number of nodes
      # indices of all columns (ravelled for 2nd and subsequent dimensions) tiled n times
      ic = np.repeat(np.arange(dim1prod, dtype=int).reshape(dim1)[None, :], n, axis=0)
      i1 = np.ravel_multi_index((i1, ic), (z.shape[0], dim1prod))
      i2 = np.ravel_multi_index((i2, ic), (z.shape[0], dim1prod))
  
      if scalar:
          i1 = np.asscalar(i1)
          i2 = np.asscalar(i2)
      else:
          # keep dimensions of the input zi
          i1 = i1.reshape(sz)
          i2 = i2.reshape(sz)
  
      return i1, i2, w1, w2
  
  # called from above
  # from pkg_interp.py
  def vinterp1d(gdepr, z, extrap_threshold):
      n = len(z)
      i2 = np.searchsorted(z, gdepr, side='right')
      ileft = i2==0  # dst layers shallower than 1st src layer
      irght = i2==n  # dst layers deeper than last src layer
      i2[ileft] = 1
      i2[irght] = n-1
      i1 = i2 - 1
      w1 = (z[i2] - gdepr) / (z[i2] - z[i1])
      w1[ileft] = 1  # this is nearest neighbour extrapolation to shallower layers
      w1[irght] = 0  # this is nearest neighbour extrapolation to deeper layers
      if extrap_threshold is not None:
          # drop points beyond extrap threshold
          invalid = np.logical_or(gdepr < z[ 0] - extrap_threshold,
                                  gdepr > z[-1] + extrap_threshold)
          w1[invalid] = np.nan
      w2 = 1 - w1
      return i1,i2,w1,w2

# class utils from etc/data_conversion_functions - GO 2022-12
class general_utils():
    def makePath(self, path):
        return os.path.abspath(os.path.expanduser(path))

    def makeCheckDirs(self, path):
        if not os.path.exists(path):
            os.makedirs(path)

        return path

    def getNextVal(self, a, pad):
        checkStr = a[a.find('.') + 1:]

        cI = 0
        while np.logical_and(cI < len(checkStr) - 1, checkStr[cI].isdigit()):
            cI = cI + 1
            if cI == 5:
                break
        cutInd = a.find('.') + 1 + cI

        startString = a[:a.find('.')]
        checkFloat = 0.
        for i in range(len(startString)):
            if startString[i].isdigit():
                continue
            else:
                startInd = i
                checkFloat = float(a[startInd:a.find('.') + 1 + cI])

        if checkFloat == 0.:
            checkFloat = float(a[:cutInd])

        if checkFloat == float(pad):
            return np.nan, cutInd
        else:
            return checkFloat, cutInd

    def makeInt(self, strL, inds):
        S = ''
        for i in inds:
            S = S + str(strL[i])

        return int(S)

    # helpers to filter profiles with possible null values -GO
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
    
    # removes nans from all if one has nan
    # assumes same shape of all arrays
    def remove_nan(T, S, P, badvals):
        ikeep = ~ (np.isnan(T)|np.isnan(S)|np.isnan(P))
        badvals += np.sum(~ikeep)
        T,S,P=T[ikeep],S[ikeep],P[ikeep]
        return T,S,P,badvals
    



def seconds_in_one_month(year, month):
    _, num_days = calendar.monthrange(year, month)
    return num_days * 24 * 60 * 60


# generic mvg avg filter using user-defined months and / or days
def apply_mvgavg_filter(time_sorted, var_sorted, mos=1, dys=30, use_nanmean=True):
    three_months_in_seconds = mos * dys * 24 * 60 * 60

    result_time = []
    result_var = []

    start_time = datetime.fromtimestamp(time_sorted[0])
    end_time = datetime.fromtimestamp(time_sorted[-1])

    current_time = start_time
    while current_time <= end_time:
        # Calculate the end of the 4-month period
        end_of_period = current_time + timedelta(seconds=three_months_in_seconds)

        # Create a time mask for the current 4-month period
        time_mask = (time_sorted >= current_time.timestamp()) & (time_sorted < end_of_period.timestamp())

        # Calculate the average var concentration for the current 4-month period
        if use_nanmean:
            average_var = np.nanmean(var_sorted[time_mask])
        else:
            average_var = np.mean(var_sorted[time_mask])

        result_time.append(current_time.timestamp())
        result_var.append(average_var)

        current_time += timedelta(seconds=three_months_in_seconds)

    # Convert the result_time back to datetime objects
    result_time = [datetime.fromtimestamp(ts) for ts in result_time]

    # Convert the result arrays to numpy arrays
    result_time = np.array(result_time)
    result_var = np.array(result_var)
    return result_time, result_var


# get stats for taylor fig from SST scores from buoys
# bu - buoy
# mv_avg - in days
# SSTbuoy_scores - from pyap (loaded pickle file output from analyse.py)
# run_sname - alt name for model run
# use_nanmean - not used!
def get_taylor_stats_SST(stn, mv_avg, SST_scores, run_sname, use_nanmean):
    # if 1 day then just use the precalculated pyap stats

    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        ccoef = scores['pearson']
        wss = scores['skill1981']

        # for target
        rmsd = scores['rmse']
        bias = scores['bias']
        nrmsd = rmsd / stdev_obs
    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o['temperature'].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m['temperature'].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        stats = sm.taylor_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, 'data')
        sdev = stats['sdev'][1]
        crmsd = stats['crmsd'][1]
        ccoef = stats['ccoef'][1]
        wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        mod_norm_stdev = stdev_mod / stdev_obs

        ncrmsd = crmsd / stdev_obs

    return [stn, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss]


def get_taylor_stats_LH(stn, mv_avg, SST_scores, run_sname, use_nanmean, var_ts):
    # get stats for taylor fig from scores from LH
    # created by HGo 2023-08-18 - only diff from above is
    #    there are both salin and temp data
    # stn - lighthouse
    # mv_avg - in days
    # SSTbuoy_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # use_nanmean - not used!

    # if 1 day then just use the precalculated pyap stats

    if var_ts == "T":
        var_ = 'temperature'
    elif var_ts == "S":
        var_ = 'salinity'
    else:
        print("issue with var passed to get_taylor_stats")

    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores'][var_]
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        ccoef = scores['pearson']
        wss = scores['skill1981']

    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o[var_].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m[var_].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        stats = sm.taylor_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, 'data')
        sdev = stats['sdev'][1]
        crmsd = stats['crmsd'][1]
        ccoef = stats['ccoef'][1]
        wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        mod_norm_stdev = stdev_mod / stdev_obs

        ncrmsd = crmsd / stdev_obs

    return [stn, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss]


def get_target_stats_SST(stn, mv_avg, SST_scores, run_sname, use_nanmean=True, augment_rmsd=True):
    # get stats for target fig from SST scores
    # stn - buoy
    # mv_avg - in days
    # SST_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # augment_rmse - multiply RMSD by sign of mod_stdev - obs
    # use_nanmean - not used!
    # if 1 day then just use the precalculated pyap stats
    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores']
        bias = scores['bias']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        # mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        rmsd = scores['rmse']
        nrmsd = rmsd / stdev_obs



    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o['temperature'].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m['temperature'].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        # returns (optionally normalized)
        # {'bias': bias, 'crmsd': crmsd, 'rmsd': rmsd}
        stats = sm.target_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, norm=True)
        bias = stats['bias']
        ncrmsd = stats['crmsd']
        nrmsd = stats['rmsd']

        # sdev = stats['sdev'][1]
        # crmsd = stats['crmsd'][1]
        # # ccoef = stats['ccoef'][1]
        # # wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        # # mod_norm_stdev = stdev_mod / stdev_obs
        # ncrmsd = crmsd / stdev_obs

    if augment_rmsd:
        if (stdev_mod - stdev_obs) < 0:
            ncrmsd = ncrmsd * -1

    return [stn, bias, ncrmsd, nrmsd]


def get_target_stats_LH(stn, mv_avg, LH_scores, var_ts, run_sname, use_nanmean=True, augment_rmsd=True):
    # get stats for target fig from LH scores
    # stn - buoy
    # mv_avg - in days
    # LH_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # augment_rmse - multiply RMSD by sign of mod_stdev - obs
    # use_nanmean - not used!

    if var_ts == "T":
        var_ = 'temperature'
    elif var_ts == "S":
        var_ = 'salinity'
    else:
        print("issue with var passed to get_taylor_stats")

    if mv_avg == 1:
        scores = LH_scores[run_sname][stn]['filt_interp_data']['scores'][var_]
        bias = scores['bias']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        # mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        rmsd = scores['rmse']
        nrmsd = rmsd / stdev_obs

    else:
        filt_interp_data_o = LH_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = LH_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o[var_].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m[var_].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        # returns (optionally normalized)
        # {'bias': bias, 'crmsd': crmsd, 'rmsd': rmsd}
        stats = sm.target_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, norm=True)
        bias = stats['bias']
        ncrmsd = stats['crmsd']
        nrmsd = stats['rmsd']

        # sdev = stats['sdev'][1]
        # crmsd = stats['crmsd'][1]
        # # ccoef = stats['ccoef'][1]
        # # wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        # # mod_norm_stdev = stdev_mod / stdev_obs
        # ncrmsd = crmsd / stdev_obs

    if augment_rmsd:
        if (stdev_mod - stdev_obs) < 0:
            ncrmsd = ncrmsd * -1

    return [stn, bias, ncrmsd, nrmsd]
