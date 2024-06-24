# created by G Oldford Aug 30 2023
# purpose: run trend detection analysis on anomaly time series, cell-wise



# Function to compute the slope of the linear trend
#def compute_slope(x, y):
#    regression = LinearRegression()
#    regression.fit(x, y)
#    return regression.coef_[0]
#
#def compute_slope_GO(x,y):
#    m, c = np.linalg.lstsq(x, y, rcond=None)[0]
#    return m, c
#
#def get_bootstrapped_CI(n_iter, ts, ts_time):
#    values = ts.values[~np.isnan(ts.values)]
#    time = ts_time[~np.isnan(ts.values)]
#    # Fit linear regression to the original data
#    # regression = LinearRegression()
#    # regression.fit(time.reshape(-1, 1), values)
#    A = np.vstack([time, np.ones(len(time))]).T
#    orig_slope, orig_bias = compute_slope_GO(A,values)
#
#    # Initialize an array to store the slopes
#    bootstrap_slopes = np.zeros(n_iter)
#    bootstrap_biases = np.zeros(n_iter)
#
#    # Perform bootstrapping
#    for i in range(n_iter):
#
#        # Generate a bootstrap sample by resampling with replacement
#        indices = np.random.choice(len(values), len(values), replace=True)
#        bootstrap_sample = values[indices]
#        bootstrap_time = time[indices]
#        A = np.vstack([bootstrap_time, np.ones(len(bootstrap_time))]).T
#
#        # Compute the slope of the linear trend on the bootstrap sample
#    #     bootstrap_slope = compute_slope_GO(bootstrap_time.reshape(-1, 1), bootstrap_sample)
#        bootstrap_slope, bootstrap_b = compute_slope_GO(A, bootstrap_sample)
#        bootstrap_slopes[i] = bootstrap_slope
#        bootstrap_biases[i] = bootstrap_b
#        
#    return bootstrap_slopes, bootstrap_biases, orig_slope, orig_bias

import numpy as np
from scipy.stats import linregress
import os
import xarray as xr
from scipy.stats import kendalltau # not used
import pymannkendall as mk
import netCDF4 as nc
from netCDF4 import Dataset

mod = '216'
base_dir = '/project/6006412/goldford/data_temp/extract_results_' + mod + '/' 


# this is to test whether Mann Kendall gives higher significance
# to weekly  or monthly sampling rather than seasonal avg anom?
use_avg = True 
avg_length = 'month' # can only by 'year' or 'month' right now, defaults to week (assumed to be from input)
alpha = 0.05 # p val thresh
use_abs = False # 20230922 - introduced to experiment with absolute anom trend (magnitude) for temp in spring


period_catch = 12 # change this to 3 if seasonal, 12 if annual files, 52 if weekly, etc

def mann_kendall_trend(data):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    #data_series = pd.Series(data)
    
    # Apply the Mann-Kendall test
    result = mk.original_test(data, alpha=alpha)
    
    return result.slope, result.h, result.s
    
# watch out b/c period can't be changed via the ufunc call
def mann_kendall_seasonal(data, period=period_catch):
    # You may need to preprocess the data here if it's not in the right format
    # For example, convert it to a pandas Series if it's not already
    #data_series = pd.Series(data)
    
    # Apply the Mann-Kendall test
    result = mk.seasonal_test(data, period=period, alpha=alpha)
    
    return result.slope, result.h, result.s

def trends_3D(variables):
 for variable in variables:
   for depth_group in depth_groups:
     for season in seasons:
         filename = f'weekly_anom_{variable}_{depth_group}_1980-2018_{season}.nc'
         anom_f = os.path.join(base_dir, 'anomalies/', filename)
         
         print("getting trend for " + variable + " in " + season + " and depth group " + depth_group)
         ds = xr.open_dataset(anom_f)  # Open the NetCDF file using xarray
         
         if use_abs:
           ds[variable] = abs(ds[variable])
         
         #tau_array = np.apply_along_axis(lambda x: kendalltau(range(len(x)), x)[0], axis=0, arr=ds)      
                  
           # Mann Kendall test 
           # pymannkendall
           # takes a vector list as input
#           slopes, trend_tf = xr.apply_ufunc(
#               mann_kendall_trend,
#               yearly_avgs,
#               input_core_dims=[["year"]],
#               output_core_dims=[[],[]],
#               vectorize=True,
#               dask="allowed"#,
#               #output_dtypes=[float]
#             )
         
         if use_avg:
           # Calculate yearly average for each cell
           # xarray does not currently allow mulitple conditions in groupby
           if avg_length == 'month':
             ds['time_counter'] = ds['time_counter'].dt.strftime('%Y-%m')
             avgs = ds[variable].groupby('time_counter').mean(dim='time_counter')
             coredim = 'time_counter'
             # Mann Kendall test 
             slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_seasonal,
               avgs,
               input_core_dims=[[coredim]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
               )
           elif avg_length == 'year':
             avgs = ds[variable].groupby("time_counter." + avg_length).mean(dim="time_counter")
             coredim = avg_length
             # Mann Kendall test 
             slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_trend,
               avgs,
               input_core_dims=[[coredim]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
               )
           else:
             print("to do: error handle")     
                         
         else:   
           # Mann Kendall test 
           slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_trend,
               ds[variable],
               input_core_dims=[["time_counter"]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
             )

         # Create a new dataset for writing
         trend_dataset = xr.Dataset(
          {
              "nav_lat": ds["nav_lat"],
              "nav_lon": ds["nav_lon"],
              "slope": (("y", "x"), slopes.values),
              "sig_tf":(("y", "x"), trend_tf.values),
              "s_val":(("y", "x"), s_val.values),
              },
          coords={"y": ds["y"], "x": ds["x"]},
         )

         # Set attributes for variables
         trend_dataset["nav_lat"].attrs = ds["nav_lat"].attrs
         trend_dataset["nav_lon"].attrs = ds["nav_lon"].attrs
         
         # Write the dataset to a new NetCDF file
         if use_abs: #experimental
           output_filename = f'anom_trend_{variable}_{depth_group}_{season}_abs.nc'
         else:
           output_filename = f'anom_trend_{variable}_{depth_group}_{season}.nc'
           
         if use_avg:
           output_path = os.path.join(base_dir, 'trend_output/anom_trend_from' + avg_length + '/', output_filename)
         else:
           output_path = os.path.join(base_dir, 'trend_output/anom_trend_fromweekly/', output_filename)
         trend_dataset.to_netcdf(output_path)

         ds.close()  # Close the dataset
         
# to do: the only diff from 3D is a nested loop for depth - combine
def trends_2D(variables):
  for variable in variables:
     for season in seasons:
     
        if season == 'annual':
          period_catch = 12
        else:
          period_catch = 3
     
        filename = f'weekly_anom_{variable}_1980-2018_{season}.nc'
        anom_f = os.path.join(base_dir, 'anomalies/', filename)
        
        print("getting trend for " + variable + " in " + season)
        ds = xr.open_dataset(anom_f)  # Open the NetCDF file using xarray
        if use_abs:
           ds[variable] = abs(ds[variable])
        
        if use_avg:
        
          if avg_length == 'month':
             ds['time_counter'] = ds['time_counter'].dt.strftime('%Y-%m')
             avgs = ds[variable].groupby('time_counter').mean(dim='time_counter')
             coredim = 'time_counter'
             slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_seasonal,
               avgs,
               input_core_dims=[[coredim]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
             )
          elif avg_length == 'year':
             avgs = ds[variable].groupby("time_counter." + avg_length).mean(dim="time_counter")
             coredim = avg_length
             slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_trend,
               avgs,
               input_core_dims=[[coredim]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
             )
          else:
             print("to do: error handle")     
                  
        else:
          # Mann Kendall test 
          # pymannkendall
          # takes a vector list as input
          slopes, trend_tf, s_val = xr.apply_ufunc(
               mann_kendall_trend,
               ds[variable],
               input_core_dims=[["time_counter"]],
               output_core_dims=[[],[]],
               vectorize=True,
               dask="allowed"#,
               #output_dtypes=[float]
             )
         
        
        #trend through each pixel-slice 
        # https://stackoverflow.com/questions/66594056/linear-regression-on-each-grid-cell-across-time-dim
        #lin_fit = yearly_avgs.polyfit('year', deg=1, skipna=True)
        # y = mx + b
        #a = lin_fit.sel(degree=1) # slopes of linregress
        #b = lin_fit.sel(degree=0)
        
        # Create a new dataset for writing
        trend_dataset = xr.Dataset(
            {
                "nav_lat": ds["nav_lat"],
                "nav_lon": ds["nav_lon"],
                "slope": (("y", "x"), slopes.values),
                "sig_tf":(("y", "x"), trend_tf.values),
                "s_val":(("y", "x"), s_val.values),
            },
            coords={"y": ds["y"], "x": ds["x"]},
        )

        # Set attributes for variables
        trend_dataset["nav_lat"].attrs = ds["nav_lat"].attrs
        trend_dataset["nav_lon"].attrs = ds["nav_lon"].attrs

        # Write the dataset to a new NetCDF file 
        if use_abs: 
          output_filename = f'anom_trend_{variable}_{season}_abs.nc'
        else:
          output_filename = f'anom_trend_{variable}_{season}.nc'
        if use_avg:
          output_path = os.path.join(base_dir, 'trend_output/anom_trend_from' + avg_length + '/', output_filename)
        else:
          output_path = os.path.join(base_dir, 'trend_output/anom_trend_fromweekly/', output_filename)
        trend_dataset.to_netcdf(output_path)

        ds.close()  # Close the dataset



years = np.arange(1980, 2019, 1)
#seasons = ['winter', 'spring', 'summer', 'fall', 'annual']
seasons = ['annual']
depth_groups = ["0to30m", "30to150m", "gt150m", "allz"]
#variables = ['votemper', 'vosaline', 'vomecrty', 'vozocrtx']
variables = ['vosaline', 'votemper']

if __name__ == "__main__":
  #seasons = ['annual']
  #depth_groups = ["allz"]
  #variables = ['vozocrtx']
  trends_3D(variables)
  #variables = ['mldkz5', 'mldr10_1']
  #trends_2D(variables)


