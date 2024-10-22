# Created by G Oldford Oct 2024
# last edit: Oct 2024
0
# Purpose: On server (here for testing), process monthly fields (outputs from ksh) in 3D to monthly depth int pacea
# Inputs: NEMO daily mean files on server pre-processed to monthly files by ksh script
#
# Output: single NC file per var per stat for period (eg 1980to2019_temperature_mean.nc)
#
# Process:
# 1- get monthly files pre-processed as means, min, max etc from NEMO daily outputs
# 2- compute the monthly means (all depths), min, max, std
# 3- clip to remove land and cells greater than max depth
# 4- clip the results to remove the Fraser River (artificial jaunt)
# 5- do depth integration or selection of depth layer
# 6- write out to NC file containing monthly data
#
# Notes / To-dos:
# 2024-10-11 - added clip of Fraser R.
# 2024-10-11 - fixes to indexing of layers

import sys
sys.path.append("C://Users//Greig//Documents//github//HOTSSea_v1_NEMOandSupportCode//desktop//code//manuscript_figs")
# noinspection PyUnresolvedReferences
from GO_tools import buildSortableString, read_sdomains # last to from GO_helpers
import netCDF4 as nc
import os
import numpy as np
from matplotlib.path import Path as Path


# ////////////////////// PATHS ////////////////////////
# local test p's
nc_p = "C:/Users/Greig/Downloads/NEMO_out/" # double fwd slashes don't work with os path join
# eg SalishSea1500-RUN216_MonthlyMean_grid_T_2D_y1980m10.nc
project_p = ""
tmask_p = "../../data/mesh mask/"
tmask_f = "mesh_mask_20210406.nc"
bathy_p = "..//..//data//bathymetry//bathy_salishsea_1500m_20210706.nc"
out_p = ""

frmask_p = "..//..//data//"
frmask_f = "FraserRiverMaskPnts.yml"
frmask_fullp = os.path.join(frmask_p, frmask_f)

# server paths
# nc_p = "DATA/SS1500-RUN216/NEMO_monthly_NC/"
# project_p = "/project/6006412/goldford/ECOSPACE/"
# tmask_p = "..//..//data//mesh mask//mesh_mask_20210406.nc"
# out_p = project_p + "DATA/SS1500-RUN{modelrun}/NEMO_monthly_NC_pacea/"

modelrun = "216"
variables = {'votemper':"temperature"}#, 'votemper','vosaline',
# for reference, the centre of the cell depths (gdept_0), not accounting for ssh stretch:
#     [[0.5000003   1.5000031   2.5000114   3.5000305   4.5000706   5.5001507
#       6.5003104   7.500623    8.501236    9.502433   10.5047655  11.509312
#       12.518167   13.535412   14.568982   15.634288   16.761173   18.007135
#       19.481785   21.389978   24.100256   28.229916   34.685757   44.517723
#       58.484333   76.58559    98.06296   121.866516  147.08946   173.11449
#       199.57304   226.2603    253.06664   279.93454   298.58588   308.9961
#       360.67453   387.6032    414.5341    441.4661]]
# and the max and min delineation of each vertical cell (deptw_0), not accounting for ssh stretch:
# [[  0.          1.0000012   2.0000064   3.0000193   4.0000467   5.000104
#     6.000217    7.0004406   8.000879    9.001736   10.003407   11.006662
#    12.013008   13.025366   14.049429   15.096255   16.187304   17.364035
#    18.705973   20.363474   22.613064   25.937412   31.101034   39.11886
#    50.963238   67.05207    86.96747   109.73707   134.34593   160.02956
#   186.30528   212.89656   239.65305   266.4952    293.3816    303.79166
#   347.2116    374.1385    401.06845   428.       ]]
depth_groups = {"surface":{"min_t_idx":0, "max_t_idx":3, "min_z":0, "max_z":4},
                "avg0to30m":{"min_t_idx":0, "max_t_idx": 21, "min_z":0, "max_z":30},
                "avg30to150m":{"min_t_idx":21, "max_t_idx": 28, "min_z":30, "max_z":150},
                "avg150mtoBot":{"min_t_idx":28, "max_t_idx": 39, "min_z":150, "max_z":"Bot"},
                "bottom":{"min_t_idx":"na", "max_t_idx":"na", "min_z":"na", "max_z":"na"} # min max ignored
                }
#depth_groups = {"avg150mtoBot":{"min_t_idx":28, "max_t_idx": 39, "min_z":150, "max_z":"Bot"}}
stat_types = ["Mean", "Min", "Max", "Std"]
stat_types = ["Min", "Max"]

start_yr = 1980
end_yr = 1980
start_mo = 1
end_mo = 12

# copy NC structure from a given input file
def create_nc_file_from_template(template_file, output_filename, var, var_name, time_len, dep_bin, stat):
    with nc.Dataset(template_file, 'r') as template:
        dataset = nc.Dataset(output_filename, 'w', format='NETCDF4')
        if dep_bin == "bottom" or dep_bin == "surface":
            dataset.description = f"Contains monthly {stat.lower()} of daily mean {var_name} from {start_yr} to {end_yr} for {dep_bin}."
        else:
            dataset.description = f"Contains monthly {stat.lower()} of daily mean {var_name} from {start_yr} to {end_yr}, depth integrated means for {dep_bin}."

        # copy dim (time is unlimited)
        for name, dimension in template.dimensions.items():
            if name == 'time_counter':
                dataset.createDimension(name, time_len)
            else:
                dataset.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

        # copy vars
        for name, variable in template.variables.items():
            if name == var:
                #fill_val = 0 # A Edwards request - not working
                fill_val = variable._FillValue
                new_var = dataset.createVariable(f'{var}', variable.datatype, ('time_counter', 'y', 'x'),
                                                 fill_value=fill_val)
                new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs() if k != "cell_methods"})
                new_var.units = variable.units

                if dep_bin == "bottom" or dep_bin == "surface":
                    if dep_bin == "bottom":
                        new_var.long_name = f"Monthly {stat.lower()} of daily mean {variable.long_name} for {dep_bin} layer."
                    else:
                        new_var.long_name = f"Monthly {stat.lower()} of daily mean {variable.long_name} for {dep_bin} (0 to 4m)."
                else:
                    new_var.long_name = f"Monthly {stat.lower()} of daily mean {variable.long_name} computed as depth-integrated means for {dep_bin}"
            else:
                # copy variables for the dims
                # don't copy the var of interest (update it later)
                if not name in ['time_counter', 'cell_area',
                                'nav_lat', 'nav_lon',
                                'nav_lon_bnds', 'nav_lat_bnds'
                                ]:
                    continue

                out_var = dataset.createVariable(name, variable.datatype, variable.dimensions)
                for k in variable.ncattrs():
                    out_var.setncatts({k: variable.getncattr(k)})
                    # out_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})

                # not needed, just double check
                if name not in ['votemper', 'vosaline', 'time_counter']:
                    out_var[:] = template.variables[name][:]

        return dataset, dataset.variables[f'{var}'], dataset.variables['time_counter']

# get masks, z levs, widths
with nc.Dataset(os.path.join(tmask_p, tmask_f), 'r') as mesh:
    # shapes of (1, 40, 299, 132)
    # tmask - land mask in 3D, t pnt
    # e3t0 - depth bin widths
    # gdept_0 - centres of the depth levs
    tmask = mesh.variables['tmask'][:]
    e3t0 = mesh.variables['e3t_0'][:]
    gdept_0 = mesh.variables['gdept_0'][:]
    gdepw_0 = mesh.variables['gdepw_0'][:]
    navlat = mesh.variables['nav_lat'][:]
    navlon = mesh.variables['nav_lon'][:]

# unsure why tmask is exporting zeros instead of nodata but this works
tmask_land = np.ma.masked_where(tmask[:,0,...] == 0, tmask[:,0,...])
tmask_fix = np.ma.masked_where(tmask == 0, tmask)


# prep array with indices of deepest lev
bottom_idxs = np.zeros([tmask.shape[2],tmask.shape[3]], dtype=np.int32)
# x_grid, y_grid = np.meshgrid(np.arange(tmask.shape[2]), np.arange(tmask.shape[3]), indexing='ij')
for i in range(tmask.shape[2]):
    for j in range(tmask.shape[3]):
        tmask2 = tmask[0, :, i, j]
        npmax = np.max(tmask2)
        if npmax == 0:
            bottom_idx = -1
        else:
            bottom_idx = np.where(tmask[0, :, i, j])[0][-1]
        bottom_idxs[i,j] = bottom_idx

# prep fraser river mask
frmask_pts = read_sdomains(frmask_fullp)
frmask = np.ones([1,tmask.shape[2],tmask.shape[3]], dtype=np.int32)
for i in range(navlat.shape[0]):
    for j in range(navlat.shape[1]):
        lat = navlat[i, j]
        lon = navlon[i, j]
        pos = [lat,lon]
        for sd in frmask_pts.keys():
            if Path(frmask_pts[sd]).contains_point(pos):
                frmask[0,i,j] = 0
frmask = np.ma.masked_where(frmask == 0, frmask)



for stat in stat_types:
    print(stat)
    template_file = os.path.join(nc_p, f'SalishSea1500-RUN216_Monthly{stat}_grid_T_y1980m01.nc')
    for var, var_name in variables.items():
        print(var, var_name)
        for dep_bin, dep_rng in depth_groups.items():
            print(dep_bin)
            min_t_idx = dep_rng['min_t_idx']
            max_t_idx = dep_rng['max_t_idx']
            min_z = dep_rng['min_z']
            max_z = dep_rng['max_z']

            # create empty output file using template
            if dep_bin == 'bottom':
                output_filename = f'hotssea_{start_yr}to{end_yr}_bottom_{var_name}_{stat.lower()}.nc'
            elif dep_bin == 'surface':
                output_filename = f'hotssea_{start_yr}to{end_yr}_surface_{var_name}_{stat.lower()}.nc'
            elif dep_bin == 'avg150mtoBot':
                output_filename = f'hotssea_{start_yr}to{end_yr}_avg{min_z}to{max_z}_{var_name}_{stat.lower()}.nc'
            else:
                output_filename = f'hotssea_{start_yr}to{end_yr}_avg{min_z}to{max_z}m_{var_name}_{stat.lower()}.nc'

            time_len = (end_yr - start_yr + 1) * (end_mo - start_mo + 1)
            nc_file, nc_var, nc_tc = create_nc_file_from_template(template_file, os.path.join(out_p, output_filename),
                                                                  var, var_name, time_len, dep_bin, stat)
            time_index = 0

            # open 3D files containing monthly means of daily output stats (mean, min, max, std) for each z
            for yr in range(start_yr,end_yr+1):
                print(yr)
                for mo in range(start_mo,end_mo+1):
                    mo = buildSortableString(mo, 2)
                    dat_f = f'SalishSea1500-RUN{modelrun}_Monthly{stat}_grid_T_y{yr}m{mo}.nc'
                    dat_p = os.path.join(nc_p, dat_f)
                    with nc.Dataset(dat_p) as dat:
                        time_cnt = dat.variables['time_counter'][:]
                        var_dat = dat.variables[var][:]
                        # navlon = dat.variables['nav_lon'][:]
                        # navlat = dat.variables['nav_lat'][:]

                    if dep_bin == "bottom":

                        dat_fin = np.empty([1,var_dat.shape[2],var_dat.shape[3]])
                        for i in range(var_dat.shape[2]):
                            for j in range(var_dat.shape[3]):
                                dat_fin[0, i, j] = var_dat[0,bottom_idxs[i,j],i,j]
                    else:
                        # depth average, applying land mask in 3D
                        dat_fin = (np.sum(tmask_fix[:,min_t_idx:max_t_idx+1,...] * e3t0[:,min_t_idx:max_t_idx+1,...] * var_dat[:,min_t_idx:max_t_idx+1,...], axis=1) /
                                   np.sum(tmask_fix[:,min_t_idx:max_t_idx+1,...] * e3t0[:,min_t_idx:max_t_idx+1,...], axis=1))

                    # no sea ice so not accurate, only occurs in a few cells and months
                    # not using <= 0 b/c it messes up fill_vals, not sure why (when viewed in R)
                    if var == "votemper":
                        dat_fin[(dat_fin > -10) & (dat_fin < 0)] = 0.01

                    # apply fraser river mask - adding tmask_land to make outputs 'nodata' instead of zeros
                    dat_fin = dat_fin * frmask * tmask_land

                    # cross checks on one cell
                    # print(var_dat.shape)
                    # print(e3t0.shape)
                    # print(var_dat[0,min_z:max_z,175,75])
                    # print(np.sum(var_dat[0,min_z:max_z,175,75] * e3t0[0,min_z:max_z,175,75]) / np.sum(e3t0[0,min_z:max_z,175,75]))
                    # print(dat_fin[0,175,75])
                    # print('dat_fin shape:')
                    # print(dat_fin.shape)

                    nc_var[time_index, :, :] = dat_fin[0, :, :]  # Assuming the shape is (1, lat, lon)

                    nc_tc[time_index] = time_cnt
                    time_index += 1

            nc_file.close()

            # crosscheck
            with nc.Dataset(output_filename) as dat:
                print(dat.description)
                # print(dat.dimensions)
                #print(dat.variables)
                navlon = dat.variables['nav_lon'][:]
                navlat = dat.variables['nav_lat'][:]
                time_cnt = dat.variables['time_counter'][:]
                # print(time_cnt)
                var_out = dat.variables[var][:]
                # print(var_out.shape)
                print(var_out[:,175,75])


