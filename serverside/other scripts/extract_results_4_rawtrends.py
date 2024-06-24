# created by G Oldford Aug 30 2023
# purpose: run trend detection analysis on raw time series, cell-wise


import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress
import os
import xarray as xr

import pymannkendall as mk

base_dir = '/project/6006412/goldford/data_temp/extract_results/' 


def mann_kendall_trend(data):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    #data_series = pd.Series(data)
    
    # Apply the Mann-Kendall test
    result = mk.original_test(data)
    
    return result.slope, result.h
    
def mann_kendall_seasonal(data, period):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    #data_series = pd.Series(data)
    
    # Apply the Mann-Kendall test
    result = mk.seasonal_test(data, period=period)
    
    return result.slope, result.h

def trend_raw_3D(variables):
 for variable in variables:
   print(variable)
   for depth_group in depth_groups:
     print(depth_group)
     for season_name, season_months in seasons.items():
       print(season_name)
       season_data_list = []
       for year in years:      
           filename = os.path.join(base_dir, 'weekly_avg/', f'weekly_avg_{variable}_{depth_group}_{year}.nc')
           ds = xr.open_dataset(filename)
           month_data = ds[variable].sel(time_counter=ds.time_counter.dt.month.isin(season_months)) 
           season_data_list.append(month_data)
           ds.close()
       combined_data = xr.concat(season_data_list, dim='time_counter')
       yearly_avgs = combined_data.resample(time_counter='1Y').mean(dim='time_counter')
       
       # Mann Kendall test - pymannkendall
       # takes a vector list as input
       print('getting slope, trend')
       slopes, trend_tf = xr.apply_ufunc(
           mann_kendall_trend,
           yearly_avgs,
           input_core_dims=[["time_counter"]],
           output_core_dims=[[],[]],
           vectorize=True,
           dask="allowed"#,
           #output_dtypes=[float]
           )

       print('writing file')
       # Create a new dataset for writing
       trend_dataset = xr.Dataset({
              "nav_lat": ds["nav_lat"],
              "nav_lon": ds["nav_lon"],
              "slope": (("y", "x"), slopes.values),
              "sig_tf":(("y", "x"), trend_tf.values),
          },
          coords={"y": ds["y"], "x": ds["x"]},
          )
         
       # Set attributes for variables
       trend_dataset["nav_lat"].attrs = ds["nav_lat"].attrs
       trend_dataset["nav_lon"].attrs = ds["nav_lon"].attrs
         
       # Write the dataset to a new NetCDF file
       output_filename = f'raw_trend_{variable}_{depth_group}_{season_name}.nc'
       output_path = os.path.join(base_dir, 'trend_output/raw_trends/', output_filename)
       trend_dataset.to_netcdf(output_path)


seasons = {
    'winter': [12, 1, 2],
    'spring': [3, 4, 5],
    'summer': [6, 7, 8],
    'fall': [9, 10, 11]
}
years = np.arange(1980, 2019, 1)
depth_groups = ["0to30m", "30to150m", "gt150m", "allz"]

if __name__ == "__main__":
  variables = ['votemper', 'vosaline', 'vomecrty', 'vozocrtx']
  trend_raw_3D(variables)
  #variables = ['mldkz5', 'mldr10_1']
  #trends_2D(variables)


