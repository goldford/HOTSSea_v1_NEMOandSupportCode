# created Aug 30 2023 by G Oldford
# purpose: generate anomalies and climatologies from weekly nemo_results

import os
import xarray as xr
import numpy as np

mod = '216'
#mod = '203'
output_dir = '/project/6006412/goldford/data_temp/extract_results_' + mod + '/' 
pattern = 'SalishSea1500-RUN' + mod + '_1d_grid_T_y{}m{:02d}.nc'


# define seasons
seasons = {
    'winter': [12, 1, 2],
    'spring': [3, 4, 5],
    'summer': [6, 7, 8],
    'fall': [9, 10, 11],
    'annual': [1,2,3,4,5,6,7,8,9,10,11,12]
}

periods = {'1980-2018': np.arange(1980, 2019, 1),
'1980-1984': [1980, 1981, 1982, 1983, 1984],
'2014-2018': [2014, 2015, 2016, 2017, 2018]
}

depth_groups = ["0to30m", "30to150m", "gt150m", "allz"]

def calculate_weekly_climatology_3D():
    climatology = {}
    
    variables = ['votemper', 'vosaline', 'vomecrty', 'vozocrtx']
    
    for period in periods:
      years = periods[period]
    
      for variable in variables:
          print("making climatology for " + variable + " for " + period)
          for depth_group in depth_groups:
              print("making climatology for " + variable + " and for " + depth_group)
              
              for season_name, season_months in seasons.items():
                print('for season ' + season_name)
                season_data_list = []
                for year in years:       
                  filename = os.path.join(output_dir, 'weekly_avg/', f'weekly_avg_{variable}_{depth_group}_{year}.nc')
                  ds = xr.open_dataset(filename)
                  month_data = ds[variable].sel(time_counter=ds.time_counter.dt.month.isin(season_months))
                  season_data_list.append(month_data)
                  ds.close()
                    
                combined_season_data = xr.concat(season_data_list, dim='time_counter')
                climatology_key = f'{variable}_{depth_group}_{period}_{season_name}'
                climatology[climatology_key] = combined_season_data.mean(dim='time_counter')
          
                # Save climatology to NetCDF
                climatology_output_filename = os.path.join(output_dir, 'climatologies/', f'weekly_climatology_{climatology_key}.nc')
                encoding = {variable: {'zlib': True, 'complevel': 4}}
                climatology[climatology_key].to_netcdf(climatology_output_filename, encoding=encoding)
                
            
def calculate_weekly_climatology_2D():
    climatology = {}
    variables = ['mldkz5', 'mldr10_1']
    
    for period in periods:
      years = periods[period]
    
      for variable in variables:
          print("making climatology for " + variable+ " for " + period)
          
          for season_name, season_months in seasons.items():
            print('for season ' + season_name)
            season_data_list = []
          
            for year in years:
              filename = os.path.join(output_dir, 'weekly_avg/', f'weekly_avg_{variable}_{year}.nc')
              ds = xr.open_dataset(filename)
              season_data_list.append(ds[variable])
              ds.close()
          
            combined_season_data = xr.concat(season_data_list, dim='time_counter')
            climatology_key = f'{variable}_{period}_{season_name}'
            climatology[climatology_key] = combined_season_data.mean(dim='time_counter')
            
            # Save climatology to NetCDF
            climatology_output_filename = os.path.join(output_dir,'climatologies/', f'weekly_climatology_{climatology_key}.nc')
            encoding = {variable: {'zlib': True, 'complevel': 4}}
            climatology[climatology_key].to_netcdf(climatology_output_filename, encoding=encoding)

  
def calculate_and_store_anomalies_3D():
  variables = ['votemper', 'vosaline', 'vomecrty', 'vozocrtx']
  for period in periods:
    years = periods[period]
    for variable in variables:
      for depth_group in depth_groups:
        print("calculating anomalies for " + variable + " and for " + depth_group)
        for season_name, season_months in seasons.items():
          print('for season ' + season_name)
          seasondepth_group_anom_data_list = [] 
          
          filename = os.path.join(output_dir, 'climatologies/', f'weekly_climatology_{variable}_{depth_group}_{period}_{season_name}.nc')
          climatology = xr.open_dataset(filename)
              
          for year in years:
            filename = os.path.join(output_dir, 'weekly_avg/', f'weekly_avg_{variable}_{depth_group}_{year}.nc')
            ds = xr.open_dataset(filename)
            weekly_data = ds[variable].sel(time_counter=ds.time_counter.dt.month.isin(season_months))
            anom_data = weekly_data - climatology[variable]
            seasondepth_group_anom_data_list.append(anom_data)
            ds.close()
              
          combined_anom_data = xr.concat(seasondepth_group_anom_data_list, dim='time_counter')
          anom_output_filename = os.path.join(output_dir,'anomalies/', f'weekly_anom_{variable}_{depth_group}_{period}_{season_name}.nc')
          encoding = {variable: {'zlib': True, 'complevel': 4}}
          combined_anom_data.to_netcdf(anom_output_filename, encoding=encoding)
          
          climatology.close()     
            
def calculate_and_store_anomalies_2D():
    variables = ['mldkz5', 'mldr10_1']
    for period in periods:
      years = periods[period]
      for variable in variables:
        anom_data_list = []
        print("calculating anomalies for " + variable)
        for season_name, season_months in seasons.items():
          season_group_anom_data_list = []
          
          filename = os.path.join(output_dir, 'climatologies/', f'weekly_climatology_{variable}_{period}_{season_name}.nc')
          climatology = xr.open_dataset(filename)
          
          print('for season ' + season_name)
          for year in years:
            filename = os.path.join(output_dir, 'weekly_avg/', f'weekly_avg_{variable}_{year}.nc')
            ds = xr.open_dataset(filename)
            weekly_data = ds[variable].sel(time_counter=ds.time_counter.dt.month.isin(season_months))
            anom_data = weekly_data - climatology[variable]
            season_group_anom_data_list.append(anom_data)
            ds.close()
              
          combined_anom_data = xr.concat(season_group_anom_data_list, dim='time_counter')
          anom_output_filename = os.path.join(output_dir,'anomalies/', f'weekly_anom_{variable}_{period}_{season_name}.nc')
          encoding = {variable: {'zlib': True, 'complevel': 4}}
          combined_anom_data.to_netcdf(anom_output_filename, encoding=encoding)
          
          climatology.close()

if __name__ == "__main__":
    climatology = calculate_weekly_climatology_3D()
    calculate_and_store_anomalies_3D()
    climatology = calculate_weekly_climatology_2D()
    calculate_and_store_anomalies_2D()

