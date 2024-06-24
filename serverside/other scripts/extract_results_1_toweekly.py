# created by G Oldford Aug 30 2023
# purpose: extract data from nemo_result and average to weekly
#          and compress files using zlib

import os
import xarray as xr
import pandas as pd
import netCDF4 as nc
import numpy as np
from multiprocessing import Pool

mod = '216'
#mod = '203'
if mod == '203':
  data_dir = '/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN' + mod + '/CDF/{}/'
else:
  data_dir = '/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN' + mod + '/CDF/' 

output_dir = '/project/6006412/goldford/data_temp/extract_results_' + mod + '/' 
pattern = 'SalishSea1500-RUN' + mod + '_1d_grid_T_y{}m{:02d}.nc'
meshm_p = '/project/6006412/mdunphy/Forcing/grid/'
meshm_f = 'mesh_mask_20210406.nc'

years = np.arange(1985,2019,1)
# non-parallelized code below
# 3D nemo output, daily averages
# average to weekly and to depth bins
depth_groups = [
    {"name": "0to30m", "min": 0, "max": 30},
    {"name": "30to150m", "min": 30, "max": 150},
    {"name": "gt150m", "min": 150, "max": None},
    {"name": "allz", "min": 0, "max": None},
]

with nc.Dataset(os.path.join(meshm_p, meshm_f)) as mesh:
      tmask=mesh.variables['tmask'][:] # 0's where depth lev exceeds
      umask=mesh.variables['umask'][:]
      vmask=mesh.variables['vmask'][:]
      fmask=mesh.variables['fmask'][:]
      e3t0=mesh.variables['e3t_0'][:] # 'widths' or 'weights' for each depth lev
      e3u0=mesh.variables['e3u_0'][:]
      e3v0=mesh.variables['e3v_0'][:]
      e3w0=mesh.variables['e3w_0'][:]

def avg_weighted_depth(dat, var, depth_var, dep_min, dep_max):
    
    #print(tmask[0,:,175,75])
    # to do: find better way
    
    #find indices closest to max and min and create mask 
    #set vars (t,u,v,w?)
    diffs_min = np.abs(dat[depth_var].values - dep_min)
    min_idx = np.argmin(diffs_min)
    dg_mask = np.zeros(dat[depth_var].shape, dtype=int)  # Initialize the mask with zeros
    
    if dep_max is None:
      dg_mask[min_idx:] = 1
    else:
      diffs_max = np.abs(dat[depth_var].values - dep_max)
      max_idx = np.argmin(diffs_max)
      dg_mask[min_idx:max_idx + 1] = 1

    
    if depth_var == 'deptht':
      e30 = e3t0
      dmask = tmask
    elif depth_var == 'depthv':
      #var_dep_trunc = dat.sel(depthv=slice(depth_group["min"], depth_group["max"]))
      e30 = e3v0
      dmask = vmask
    elif depth_var == 'depthu':
      e30 = e3u0
      dmask = umask
    elif depth_var == 'depthw':
      e30 = e3w0
      dmask = fmask
      
    dg_mask = dg_mask[np.newaxis, :, np.newaxis, np.newaxis]
    masked_dataarray = dat.copy()
    masked_dataarray = masked_dataarray.where(dmask, drop=False, other=0) # 0's where z exceeds water col max
    masked_dataarray = masked_dataarray.where(dg_mask, drop=False, other=0) # 0's where z is not within group range
    
    e30_2 = e30 * dmask * dg_mask
    e30_2_sum = np.sum(e30_2, axis=1, keepdims=True)
    e30_weights = e30_2 / e30_2_sum
    weighted_data = masked_dataarray * e30_weights # weight the values by depth bin widths

    final_avg = np.sum(weighted_data, axis=1, keepdims=True)  
    
    return final_avg


def process_year_3D(year):

    var_dict_3d = {#'votemper': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_T_y{}m{:02d}.nc', 'depth':'deptht'},
                   #'vosaline': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_T_y{}m{:02d}.nc', 'depth':'deptht'},
                   #'vomecrty': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_V_y{}m{:02d}.nc', 'depth':'depthv'},
                   #'vozocrtx': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_U_y{}m{:02d}.nc', 'depth':'depthu'},
                   'vovecrtz': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_W_y{}m{:02d}.nc', 'depth':'depthw'}
                                                             
                  }
                  
    print("working on year " + str(year))
    for variable in var_dict_3d.keys():
        data_list = []
        pattern = var_dict_3d[variable]['pattern']
        print(variable)
        depth_var = var_dict_3d[variable]['depth']
        
        for month in range(1, 13):
            if mod == '203':
              filename = os.path.join(data_dir.format(year), pattern.format(year, month))
            else:
              filename = os.path.join(data_dir, pattern.format(year, month))
            ds = xr.open_dataset(filename)
            
            # Select the variable and depth levels
            var_data = ds[variable]
            depths = ds[depth_var]
            
            ds.close()
            
            data_list.append(var_data)
        
        # Combine data from all months
        combined_data = xr.concat(data_list, dim='time_counter')
        
        # Calculate weekly averages
        weekly_avg = combined_data.resample(time_counter='1W').mean(dim='time_counter')
        weekly_avg = weekly_avg.round(2)
        
        # Save to compressed NetCDF
        output_filename = os.path.join(output_dir,'weekly_avg/' , f'weekly_avg_{variable}_{year}.nc')
        encoding = {variable: {'zlib': True, 'complevel': 4}}
        weekly_avg.to_netcdf(output_filename, encoding=encoding)
        
        # Vertical averaging by depth groups
        for depth_group in depth_groups:
            print(depth_group)
            
            avg_var = avg_weighted_depth(weekly_avg, variable, depth_var, depth_group["min"], depth_group["max"])
            
            if depth_var == 'deptht':
              avg_var = avg_var.isel(deptht = 0)
            elif depth_var == 'depthv':
              avg_var = avg_var.isel(depthv = 0)
            elif depth_var == 'depthu':
              avg_var = avg_var.isel(depthu = 0)
            elif depth_var == 'depthw':
              avg_var = avg_var.isel(depthw = 0)
            else:
              print('problem with dropping depth dim.')
            
            # write out
            depth_group_name = depth_group["name"]
            depth_group_output_filename = os.path.join(output_dir,'weekly_avg/', f'weekly_avg_{variable}_{depth_group_name}_{year}.nc')
            avg_var.to_netcdf(depth_group_output_filename, encoding={variable: {'zlib': True, 'complevel': 4}})



def process_year_2D(year):
    
    var_dict_2D = {'mldkz5': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_T_2D_y{}m{:02d}.nc'},
                   'mldr10_1': {'pattern':'SalishSea1500-RUN' + mod + '_1d_grid_T_2D_y{}m{:02d}.nc'}                                           
                  }
    
    print("working on year " + str(year))
    for variable in var_dict_2D.keys():
        data_list = []
        
        print(variable)
        pattern = var_dict_2D[variable]['pattern']
        
        for month in range(1, 13):
        
            #filename = os.path.join(data_dir.format(year), pattern.format(year, month))
            if mod == '203':
              filename = os.path.join(data_dir.format(year), pattern.format(year, month))
            else:
              filename = os.path.join(data_dir, pattern.format(year, month))
            
            ds = xr.open_dataset(filename)
            
            # Select the variable
            var_data = ds[variable]
            ds.close()
            data_list.append(var_data)
        
        # Combine data from all months
        combined_data = xr.concat(data_list, dim='time_counter')
        
        # Calculate weekly averages
        weekly_avg = combined_data.resample(time_counter='1W').mean(dim='time_counter')
        #weekly_avg = weekly_avg.round(2)
        
        # Save to compressed NetCDF
        output_filename = os.path.join(output_dir,'weekly_avg/', f'weekly_avg_{variable}_{year}.nc')
        encoding = {variable: {'zlib': True, 'complevel': 4}}
        weekly_avg.to_netcdf(output_filename, encoding=encoding)


if __name__ == "__main__":
    for year in years:
        print(year)
        process_year_3D(year)
        #process_year_2D(year)
        

# parallelized code below - using daskimport os
# didn't work - dask needs min python 3.9
#import dask
#import dask.array as da
#
#def process_file(year, month):
#    filename = os.path.join(data_dir, pattern.format(year, month))
#    ds = xr.open_dataset(filename)
#    
#    # Selecting the variable and depth levels
#    votemper = ds['votemper']
#    deptht = ds['deptht']
#    
#    ds.close()
#    
#    return votemper
#
#def process_year(year):
#    tasks = []
#
#    for month in range(1, 13):
#        tasks.append(dask.delayed(process_file)(year, month))
#    
#    # Combine data from all months using Dask
#    combined_data = da.concatenate(tasks, axis=0)
#    
#    # Calculate weekly averages using Dask
#    weekly_avg = combined_data.resample(time_counter='1W').mean(dim='time_counter')
#    
#    # Reduce precision to two decimal places using Dask
#    weekly_avg = weekly_avg.round(2)
#    
#    # Save to compressed NetCDF using Dask
#    output_filename = os.path.join(output_dir, f'weekly_avg_{year}.nc')
#    weekly_avg.to_netcdf(output_filename, compute='threads', encoding={'votemper': {'zlib': True, 'complevel': 4}})
#
#if __name__ == "__main__":
#    for year in years:
#        process_year(year)

# parallelized code below - using multiprocessing lib
# throws a cross-thread error
#def process_file(year, month):
#    filename = os.path.join(data_dir, pattern.format(year, month))
#    ds = xr.open_dataset(filename)
#    
#    # Selecting the variable and depth levels
#    votemper = ds['votemper']
#    deptht = ds['deptht']
#    
#    ds.close()
#    
#    return votemper
#
#def process_year(year):
#    pool = Pool()
#    months = range(1, 13)
#    data_list = pool.starmap(process_file, [(year, month) for month in months])
#    pool.close()
#    pool.join()
#    
#    # Combine data from all months
#    combined_data = xr.concat(data_list, dim='time_counter')
#    
#    # Calculate weekly averages
#    weekly_avg = combined_data.resample(time_counter='1W').mean(dim='time_counter')
#    
#    # Reduce precision to two decimal places
#    weekly_avg = weekly_avg.round(2)
#    
#    # Save to compressed NetCDF
#    output_filename = os.path.join(output_dir, f'weekly_avg_{year}.nc')
#    encoding = {'votemper': {'zlib': True, 'complevel': 4}}
#    weekly_avg.to_netcdf(output_filename, encoding=encoding)
#        
#def main():
#    for year in years:
#        process_year(year)
#
#if __name__ == "__main__":
#    main()