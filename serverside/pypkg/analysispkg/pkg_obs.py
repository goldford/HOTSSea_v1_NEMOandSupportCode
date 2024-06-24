import fnmatch
import gsw
import netCDF4 as nc
import numpy as np
import os

from analysispkg import pkg_data
from analysispkg import pkg_eos
from analysispkg import pkg_geo
from analysispkg import pkg_utils


# OBS ARCHIVE standard variable names:

# time        - time vector
# temperature - recorded insitu temperature
# cTemp       - conservative temperature  (TEOS-10)
# aSal        - absolute salinitu (TEOS-10)
# salinity    - practical salinity (EOS-80)
# pTemp       - potential temperature     (EOS-80)

# As of March 2022:
# - Most files omit pTemp (CTD), so we calculate it at load time below.
#   It would be better to move this conversion to the obs prep phase to keep this loading code simple.
# - The rename function should be removed, it is just adding complexity.


def load_observation(ncfile, instr_type, opt=None):
    """
    Load an observation from a nc file into a dictionary
    Parameters
    ----------
    ncfile - path to a netcdf file from observations archive
    instr_type - type of instrument (TG, CTD, etc)
    opt - options yaml, only needed for HTG

    Returns
    -------
    A dictionary with the observation data

    Also some special treatment based on instrument type below

    """
    data = {'kind': instr_type, "filename":ncfile}

    _, filename = os.path.split(ncfile)

    # Possible names for time variable
    time_vars = ['time', 'time_counter']

    # Generic NC reader reads all variables
    with nc.Dataset(ncfile) as ncf:
        
        data['station'] = get_station_label(ncf)
        # TODO Remove ad-hoc station names for individual instruments below. Check if this doesn't break the analysis

        # Loop over all variables and load the time vector, note the time dimname, and start a list of time series varnames
        for varname in ncf.variables.keys():
            if varname in time_vars:
                data['time'] = pkg_data.read_time(ncf, timename=varname)
                timedimname = ncf.variables[varname].dimensions[0]
                data['timeseries'] = ['time']
                
        # Now loop over all vars and load them
        for varname in ncf.variables.keys():
             
            # Skip time as we've read it above.
            if varname in time_vars:
                continue
            
            # added to catch the empty depths in nanoose data - GO 20230612
            if varname == 'depth':
              if not ncf.variables[varname][:].any():
                continue
            
            if 'single' in ncf.variables[varname].dimensions:
                # Single values are imported as single values
                data[varname] = ncf.variables[varname][0]
            elif any(x in ncf.variables[varname].dimensions for x in ['string', 'string1', 'string2']):
                data[varname] = ncf.variables[varname][:].tostring().decode('utf-8')
            else:
                # Multiple values imported as arrays
                data[varname] = ncf.variables[varname][:]
            # Fill masked arrays, change the fill value to np.nan
            if isinstance(data[varname], np.ma.MaskedArray):
                data[varname] = pkg_data.maskedarray_to_ndarray_with_nans(data[varname])
            
            # If this variable has a time dimension then we assume it to be a timeseries, and add it to the list of
            # time series varnames. This is so we know which variables can be time-truncated by truncate_obs_timeseries()
            if 'time' in data.keys() and data['time'].size > 1 and timedimname in ncf.variables[varname].dimensions:
                data['timeseries'] += [varname]
           
        # If available, grab the "Name" attribute for use in titles etc
        # Can in theory allow for more permutations (shortname, short_name, etc)
        # LABEL and STATION are for tide gauges
        namevars = ['name', 'Name', 'shortname', 'short_name', 'Short_Name', 'LABEL', 'STATION']
        for n in namevars:
            if n in ncf.ncattrs():
                data['shortname'] = ncf.getncattr(n)
                break 
            
    # Record start and end times
    if 'time' in data.keys():
        data['tstart'] = data['time'].min()
        data['tend'] = data['time'].max()

    # Rename function to change varnames in the obs file to different names in the obs dicts.
    # This isn't ideal but if we do away with this then we have to revisit a lot of functions in the main package
    # TODO; medium priority
    def rename(data, old, new):
        # change varname; straightforward
        data[new] = data.pop(old)
        # also change it in the timeseries list if needed
        if 'timeseries' in data.keys() and old in data['timeseries']:
            data['timeseries'] = [x for x in data['timeseries'] if x != old]
            data['timeseries'] += [new]

    # TODO: remove this
    if 'longitude' in data.keys():
        rename(data, 'longitude', 'lon')
    if 'latitude' in data.keys():
        rename(data, 'latitude', 'lat')

    if instr_type == 'SST':
        data['station'] = filename.split("_SST")[0]
        data['z'] = 0

    if instr_type == 'LH':
        data['station'] = filename.split("_LH")[0]
        data['z'] = 0

    if instr_type == "CTD":
    
        #==== MK: QC to catch erroneous data (e.g. bogus 195 dbar values in Hakai CTDs)
        def reindex_data(idx):
            # Helper to extract valid values or sort
            data['Pres'] = data['Pres'][idx]
            data['temperature'] = data['temperature'][idx]
            data['cTemp'] = data['cTemp'][idx]
            data['salinity'] = data['salinity'][idx]
            data['aSal'] = data['aSal'][idx]
            data['pDen'] = data['pDen'][idx]

        # make sure values are within some valid range
        ival = (0 <= data['salinity']) & (data['salinity'] < 40) & \
               (0 <= data['temperature']) & (data['temperature'] < 40)
        if np.any(~ival):
            reindex_data(ival)

        # make sure the profile is sorted by depth
        if np.any(np.diff(data['Pres']) < 0):
            isort = np.argsort(data['Pres'])
            reindex_data(isort)
        #==== MK: end QC
        
        data['station'] = filename.split(".nc")[0]
        data['z'] = -1 * gsw.z_from_p(data['Pres'], data['lat'])

        # TODO: this should ultimately be moved to data prep phase.
        data['pTemp'] = gsw.pt_from_t(data['aSal'], data['temperature'], data['Pres'], p_ref=0.)
        data['timeseries'] += ['pTemp']

    return data


def get_station_label(ncf):
    """ Generate station ID from observation file
    ncf : netCDF file object
    """
    if "LABEL" in ncf.ncattrs():
        id_code = ncf.getncattr("LABEL")
    elif "stn_id" in ncf.ncattrs():  # HTG
        id_code = ncf.getncattr("stn_id")
    else:
        # Synthesize a label from filename
        f1 = os.path.basename(ncf.filepath())
        f2, e2 = os.path.splitext(f1)
        id_code = f2.replace(" ", "_")
    return id_code


def truncate_obs_timeseries(obs,tstart,tend):
    if 'timeseries' in obs.keys():
        c1 = np.where(np.logical_and(obs['time'] >= tstart, obs['time'] <= tend))[0]
        if len(c1) == 0:
            return None
        else:
            for varname in obs['timeseries']:
                obs[varname] = obs[varname][c1,...]
    obs['tstart'] = obs['time'].min()
    obs['tend'] = obs['time'].max()

    return obs

def truncate_mod_timeseries(mod,tstart,tend, varlist, exclude_end=False):
    if mod is None:
        return None
    if exclude_end:
        c1 = np.logical_and(mod['time'] >= tstart, mod['time'] < tend)
    else:
        # This is slow since it requires a large number of datetime object comparisons:
        #   c1 = np.logical_and(mod['time'] >= tstart, mod['time'] <= tend)
        # But since time is sorted we can use a bisection based search:
        c1 = np.zeros(mod['time'].shape, dtype=bool)
        i1 = np.searchsorted(mod['time'],tstart, side='left' )
        i2 = np.searchsorted(mod['time'],tend  , side='right')
        c1[i1:i2] = True
    if not np.any(c1):
        return None
    else:
        for varname in varlist:
            mod[varname] = mod[varname][c1,...]
    mod['tstart'] = mod['time'].min()
    mod['tend'] = mod['time'].max()

    return mod


def observations_find_nearest_point_and_waterdepth(stations, coordfile, bathyfile):
    """
    Find the indices for the nearest T-point and water level at that location

    Input is list of dict, output is list-of-dict with new fields (ilon,ilat,wd)
    """

    # Coordinate file
    with nc.Dataset(coordfile, "r") as ncid:
        glon = ncid.variables["nav_lon"][:, :]
        glat = ncid.variables["nav_lat"][:, :]
        glamf = np.squeeze(ncid.variables["glamf"][..., :, :])
        gphif = np.squeeze(ncid.variables["gphif"][..., :, :])
    glamfe, gphife = pkg_geo.expandf(glamf, gphif)

    # Bathymetry file
    with nc.Dataset(bathyfile, "r") as ncid:
        bathy = ncid.variables["Bathymetry"][:, :]
    lm = np.where(bathy > 0, 1, 0)

    out_stations = []
    for stn in stations:
        lon, lat = stn["lon"],stn["lat"]

        ilon, ilat, dist = pkg_data.find_nearest_point_fast(
            lon, lat, glon, glat, lm, glamfe, gphife
        )

        if np.isnan(ilat):
            continue
        else:
            # Get the water depth at that point
            wd = bathy[ilat, ilon]
            stn["ilon"] = int(ilon)
            stn["ilat"] = int(ilat)
            stn["wd"] = float(wd)
        out_stations += [stn]

    return out_stations


def observations_filter_time(stations, mstart, mend):
    """

    Parameters
    ----------
    stations - list of stations
    mstart - start time
    mend - end time

    Returns
    -------
    Removes stations that do not overlap the time span [mstart,mend]
    """
    out_stations = [stn for stn in stations if
                    (stn['time'][-1] >= mstart) & (stn['time'][0] <= mend)]
    return out_stations


def station_locs(stations):
    # Build list of station locations
    lons = []
    lats = []
    idx = []  # station index
    for k, stn in enumerate(stations):
        if np.isscalar(stn["lon"]):  # scalar coords
            lons += [stn["lon"]]
            lats += [stn["lat"]]
            idx += [k]
        else:  # list of coords
            lons += stn["lon"]
            lats += stn["lat"]
            idx += [k] * len(stn["lon"])
    return np.array(lons), np.array(lats), np.array(idx)


def observations_filter_location(stations, coordfile):
    """ Filters the observations to exclude stations outside of the domain

    This function is suitable both for single-position measurements (TG, CTD, ...) and
    for multi-point measurements. For multi-point measurements,
    keeps measurements with ALL points inside the domain.
    
    Input is list of dict, output is list-of-dict
    """

    # Prepare polygon following domain boundaries
    bbox = pkg_geo.domain_polygon(coordfile)

    # Build list of station locations
    lons, lats, idx = station_locs(stations)
    # Keep stations inside the domain
    flags = bbox.contains_points(np.c_[lons, lats])
    exclude = np.unique(idx[~flags])  # exclude stations with ANY point outside the domain
    include = set(range(len(stations))) - set(exclude)  # include stations that we do not exclude
    out_stations = [stations[k] for k in include]

    # The code above is equivalent to the list-comprehension below, but it is
    # faster because we call contains_points once with all of the points
    # instead of n times for each station
    # out_stations = [stn for stn in stations if np.all(bbox.contains_points(np.c_[stn["lon"], stn["lat"]]))]

    return out_stations


def observations_filter_exclude(stations, excludes):
    """ Applies list of exclude patterns and removes matches

    Input is list of dict, output is list-of-dict
    """

    def keep_station(stn, patterns):
        for pattern in patterns:
            if fnmatch.fnmatch(stn['filename'], pattern):
                #print("dropping ", stn['filename'])
                return False
        #print("keeping ", stn['filename'])
        return True

    out_stations = [stn for stn in stations if keep_station(stn, excludes)]
    return out_stations


def observations_load_and_filter(opt, instr_type, ana_start, ana_end, filter_location=True):
    """ Load observations index and filter by location, time, and exclude list.

    Parameters
    ----------
    opt
    instr_type
    ana_start
    ana_end
    filter_location : bool
        Currently set to False in CODAR and FERRY processing

    Returns
    -------

    """
    observations = pkg_utils.load_index(opt['obs'], instr_type)
    print("{} observations found".format(len(observations)))
    if filter_location:
        observations = observations_filter_location(observations, opt['file_coord'])
        observations = observations_filter_polygon(opt, observations)
    # n1 = len(observations)
    if instr_type != 'HTG':
        observations = observations_filter_time(observations, ana_start, ana_end)
    if 'exclude' in opt['extract'][instr_type].keys():
        observations = observations_filter_exclude(observations, opt['extract'][instr_type]['exclude'])
    # n2 = len(observations)
    # print("Using {} of {} {} records for current period".format(n2, n1, instr_type))
    print("{} observations to process after filtering".format(len(observations)))
    if opt.get('test', False):
        nt = 6
        if len(observations) > nt:
            observations = observations[:nt]
        print("{} observations to process (TEST MODE)".format(len(observations)))
    return observations
