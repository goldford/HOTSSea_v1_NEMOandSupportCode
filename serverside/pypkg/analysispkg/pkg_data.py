import csv
import datetime
import fnmatch
import glob
import netCDF4
import numpy as np
import os
import pickle
import scipy.io
from collections import Counter

from analysispkg import pkg_geo
from analysispkg import pkg_scan
from analysispkg import pkg_utils



def list_source_files(source_output_dir, source_output_pattern):
    """ List available model output files.
    """
    # list files
    search_str = os.path.join(source_output_dir, '**', source_output_pattern)
    fl = glob.glob(search_str, recursive=True)
    if len(fl) == 0:
        raise RuntimeError('No data files found. Search path: \n{}'.format(search_str))
    return fl


def list_source_files_sorted(opt, pattern, timename, return_times=False):
    """ List available model output files. Sorted by increasing time using first time record to sort on.
    """
    # get scanned time data
    timedata = pkg_scan.load_netcdf_index(opt)

    # get files matching the pattern
    files = fnmatch.filter(timedata.keys(),pattern)

    # get times for each file
    times = [pkg_scan.get_time(timedata, file, timename) for file in files]

    # get first time record from each file
    # t0 = [pkg_scan.get_time(timedata, file, timename)[0] for file in files]
    t0 = [timesf[0] for timesf in times]

    # sort by time
    isort = np.argsort(t0)
    files = np.array(files)[isort].tolist()
    times = [times[i] for i in isort]

    if return_times:
        return files, times
    else:
        return files


def list_source_files_sorted_in_range(opt, pattern, timename, t_start, t_end):
    """ List available model output files that have times between t_start and t_end.
    """
    files0 = list_source_files_sorted(opt, pattern, timename)

    # get scanned time data
    timedata = pkg_scan.load_netcdf_index(opt)

    def inrange(file, t0, t1):
        t = pkg_scan.get_time(timedata, file, timename)
        check = np.logical_and(t >= t0, t <= t1)
        return np.any(check)

    files = [file for file in files0 if inrange(file, t_start, t_end)]

    return files

def nemo_files_and_times(opt, pattern, timename):
    """ Time series for NEMO output with file names and time indices within each file.

    Picks unique time stamps using first occurrence in the sorted time series.
    TODO: we should raise an error when duplicates found

    Parameters
    ----------

    opt: dict
        Options dict
    pattern : str
        Glob string for file search.
    timename : str
        time vector name

    Returns
    -------

    times : list of datetime
        Output times.
    files : list of str
        Name of the file for each time, len(times).
    kt : list of int
        For each time: Index of this time step in the corresponding file; len(times).
    """
    # get scanned time data
    timedata = pkg_scan.load_netcdf_index(opt)

    # get files matching the pattern
    files = fnmatch.filter(timedata.keys(),pattern)

    mtimes, mfiles, mkt = [], [], []
    for file in files:
        t = pkg_scan.get_time(timedata, file, timename)
        nt = len(t)
        # extend global lists
        mtimes.extend(t)
        mfiles.extend([file] * nt)  # replicate file name nt times and append to global list
        mkt.extend(range(nt))  # time indices

    # exclude duplicates: np.unique chooses the values which occur first in list
    mtimes, ui = np.unique(mtimes, return_index=True)  # returns list sorted by time
    # select unique times
    mfiles = [mfiles[i] for i in ui]
    mkt = [mkt[i] for i in ui]
    return mtimes, mfiles, mkt


def read_time(ncid, timename='time', calendar='standard', offset=0):
    """Reads the time vector from a netcdf file and converts it to a datetime object.

    Used to read both obs and model files.

    Arguments:
    ncid: An opened netcdf file handle
    timename: name of the time vector; default "time"
    calendar: If calendar is not included in the netcdf attributes, use this value (default standard)
    offset: offset in seconds to add to the time vector.
    """
    if timename not in ncid.variables:
        raise RuntimeError('Time variable {} not found.'.format(timename))

    time_raw = ncid.variables[timename][:].filled()
    units = ncid.variables[timename].units
    try:
        cal = ncid.variables[timename].calendar
    except:
        cal = calendar

    if 'seconds since' in units:
        # Convert reference time to netcdf4's datetime, then to datetime.datetime
        time00 = netCDF4.num2date(0, units, cal)
        time0 = datetime.datetime(time00.year, time00.month, time00.day, time00.hour, time00.minute, time00.second)
        # Create ufunc for converting seconds to dt.timedelta
        to_timedelta = np.frompyfunc(lambda s: datetime.timedelta(seconds=s),1,1)
        # Apply ufunc to time_raw and add reference time to get array of dt.datetime
        time_var = to_timedelta(time_raw) + time0

    else:  # old method
        # convert the numbers from the file to datetime objects
        time_var = netCDF4.num2date(time_raw, units, cal)
        # Convert to array of datetime.datetime
        time_var = np.array([datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in time_var])
    # Add offset
    time_var = apply_time_offset_seconds(time_var, offset)
    return time_var


def apply_time_offset_seconds(dates, offset):
    """
    Adds offset in seconds to each datetime entry in dates
    """
    if offset == 0:
        return dates
    dt = datetime.timedelta(seconds=offset)
    tmp = [x + dt for x in dates]
    if isinstance(dates,list):
        return tmp
    elif isinstance(dates, np.ndarray):
        return np.array(tmp)
    else:
        print("apply_time_offset_seconds: unable to process dates")



def obs_to_picklefile(opt, obsfilename):
    """ Generate path and name for extracted model data file from the corresponding obs file.

    Parameters
    ----------
    opt : dict
        Options.
    obsfilename : str
        Observations data file name.

    Returns
    -------
    pickledir : str
        Where to store extracted model data.
    picklefile : str
        Name of the file with the extracted model data.

    """
    obsroot = opt['obs']['root']
    extractroot = opt['dir_extract']
    obsrelpath = os.path.relpath(obsfilename, obsroot)
    picklefile = os.path.join(extractroot, obsrelpath).replace(".nc",".pickle")
    pickledir, _ = os.path.split(picklefile)
    return pickledir, picklefile


def load_extraction(opt,obs):
    """
    Saves a model extraction 'mod' corresponding to observation 'obs' from a pickle file
    """
    pdir, pfile = obs_to_picklefile(opt, obs['filename'])
    if os.path.exists(pfile):
        with open(pfile, 'rb') as fid:
            mod = pickle.load(fid)
        return mod, pfile
    else:
        return None, pfile

def save_extraction(opt,obs,mod):
    """
    Saves a model extraction 'mod' corresponding to observation 'obs' in a pickle file

    Parameters
    ----------
    opt
    obs
    mod

    Returns
    -------
    Pickle file name
    """
    pdir, pfile = obs_to_picklefile(opt, obs['filename'])
    os.makedirs(pdir, exist_ok=True)
    mod['filename'] = pfile
    with open(pfile, 'wb') as fid:
        pickle.dump(mod, fid)
    return pfile


def parse_lcn (loc, n=2):
    """Utility: Takes a (masked) location tuple and returns a string with the values rounded to n decimal places."""
    
    lon = loc[0]
    lat = loc[1]
    lon = round(float(lon),n)
    lat = round(float(lat),n)
    loc_str = str(lat) + ' N, ' + str(abs(lon)) + ' W'

    return loc_str



def read_csv (filename, delim=','):
    """Utility: reads a csv file and returns a list of stations"""

    list_stations = []

    with open(filename, newline='') as csvfile:
        rdr = csv.reader(csvfile, delimiter=delim, quotechar='|')
        for r in rdr:
            list_stations.append([q.strip() for q in r])

    return list_stations



def trim (var, vartime, startdate, enddate):
    """Utility:  LIKELY BROKEN - DO NOT USE FOR NOW"""

    #sd=time.mktime(datetime.datetime.strptime(startdate, "%Y%m%d").timetuple())
    #print (startdate, sd, vartime[0])

    sd = pkg_utils.numdate (startdate)
    ed = pkg_utils.numdate (enddate)    

    # Grab only the data  in the window.  Everything is already in UTC.
    ind = np.where ((vartime >= sd) & (vartime < ed))
    var = var[ind]
    vartime = vartime[ind]

    return var,vartime


def nan_threshold (var, threshold):
    """Utility: Sets all values in var > threshold and < -threshold to np.nan"""

    ind = np.where (((var > threshold) | (var < -threshold)) & (np.isnan(var) == False)  )
    var[ind] = np.NaN

    return var

def replace_nans (var):
    """Utility: replaces nans in a variable with the nan-ignoring mean of the variables; also returns locations of nans."""

    nm = np.nanmean(var)
    inds = np.where (np.isnan(var))
    var[inds] = nm

    return var, nm, inds

def remove_mean (var):
    """Utility: Remove the mean from the variable. NaNs are ignored."""

    mm = np.nanmean(var)
    var[:] = var[:] - mm

    return var, mm

def missing (var, ts=(0,0,0)):
    """Utility: Calculates the percentage of a variable that is np.nan"""

    # If start,end,dt given, calculate how long the time series should be and use that as the length
    if ts != (0,0,0):
        start, end, dt = ts
        nt = (end - start).total_seconds() / dt
    else:
        nt = len(var)
    pm = np.isnan(var).sum() / nt
    pm = round(100.0*pm,3)

    return pm


def regularize (var, vartime, regtime):
    """Utility: Takes a possibly gappy or irregular time array as well as an array with consistent spacing and returns the intersection of the two."""

    regvar = np.NaN * np.zeros_like(regtime, dtype='float')

    # Ignore microseconds.  
    #In general, this should probably be "round to nearest minute"
    vartime = np.asarray([ datetime.datetime.strptime( t.strftime('%Y%m%d %H:%M:%S'), '%Y%m%d %H:%M:%S') for t in vartime] )

    _,inda, indb = np.intersect1d (vartime, regtime, return_indices=True)
    regvar[indb] = var[inda]

    return regvar


def find_nearest_point (lon,lat,glon,glat,mask,dx):
    """Utility: Finds the indicies of a station.
                lon and lat are the coordinates of the station
                glon and glat are the fields of coordinates from coord file
                    mask is the land mask, since if the closest point is on land it is not valid
                dz is the radius in degrees around which to search"""

    # This funciton only supports lon/lat values that are single floats, not arrays
    # For multiple points, call this function repeatedly (once per point).  

    #lon=lon-0.04   # no clue why I put this in at some point

    # Find the points within dx degrees of the specified lon,lat
    b=np.nonzero( (glon[:,:] < lon+dx) & (glon[:,:] > lon-dx) & (glat[:,:] < lat+dx) & (glat[:,:]> lat-dx) & (mask[:,:]>0))
    if (len(b[0]) ==0) :
        return (np.nan, np.nan, np.nan)

    npts = b[0].shape[0]    # how many points are in the range?
    dist = np.zeros(npts)   # initialize variable

    for n in range (npts):
        # Get the indices in b in a shorter form
        # note that since this uses netcdf 3 files and the scipy.io import regime, dimensions go (y,x) rather than (x,y).
        # Make sure this is consistent with the calling program - don't use without verifying that this is sensible
        # If you change the order here, change it a few lines down when you're getting the minimum distance too.
        ix = b[1][n];  iy = b[0][n]

        dist[n] = haversine(lon, lat, glon[iy,ix], glat[iy,ix])

    # Now get the minimum distance.
    ind= np.argmin(dist)            # get the index of the minimum value
    ix = b[1][ind]; iy = b[0][ind]        # get the indicies of the relevant point in the search regions

    # Return a named tuple with all the info: output

    #nearest_pt = collections.namedtuple('nearest_pt', ['id','grp', 'glon', 'glat', 'mask', 'dist', 'ix', 'iy'])
    p = (ix,iy, dist[ind]) #nearest_pt (id,grp,glon[iy,ix], glat[iy,ix], mask[iy,ix], dist[ind], ix, iy)
    return p


def haversine(lon1, lat1, lon2, lat2, rearth=6378.137):
    """ This is borrowed from from salishsea_tools
        Distance returned is in kilometres
        Operands can be arrays, but should be broadcastable.
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = rearth * c  # from m_lldist
    #km =# 6367 * c # from ss tools
    #km = 6371 * c  # from earthdist.m
    return km


def find_nearest_point_fast(lon,lat,model_lons,model_lats,mask, glamfe, gphife, max_boxes=6):

    # Find grid box containing the target point
    ie,je = pkg_geo.search(glamfe,gphife,lon,lat)

    # Define search box to find the nearest /water/ point
    i1 = max(0, ie - max_boxes)
    i2 = min(ie + max_boxes, model_lons.shape[1])
    j1 = max(0, je - max_boxes)
    j2 = min(je + max_boxes, model_lons.shape[0])

    # Prepare lists of test points in search box and their indices
    lons = model_lons[j1:j2,i1:i2].flatten()
    lats = model_lats[j1:j2,i1:i2].flatten()
    masks = mask[j1:j2,i1:i2].flatten()
    i_list = np.array([x for x in range(i1,i2)]*(j2-j1))
    j_list = j1 + np.array([x for x in range(0,(j2-j1)*(i2-i1))]) // (i2-i1)

    # Test all of these points, return closest point
    dists = haversine(
        np.array([lon] * i_list.size), np.array([lat] * j_list.size),
        lons, lats)
    dists += (1-masks)*1e20  # ensure land values fail to be closest point
    n = dists.argmin()
    if dists[n] >= 1e20:
        return np.nan, np.nan, np.nan
    j, i = j_list.item(n), i_list.item(n)
    return i, j, dists[n]


def maskedarray_to_ndarray_with_nans(fld):
    """
    Converts a np.ma.MaskedArray to a regular numpy ndarray and fills masked values with np.nan
    """
    if isinstance(fld, np.ma.MaskedArray):
        
        fld_filled = fld.filled()
        
        if fld.mask is not np.bool_(False) and fld.mask.size == fld_filled.size:
            fld_filled[fld.mask] = np.nan
        
        return fld_filled
    else:
        return fld


# code from ttide_py
# Reimplementation of the MatPlotLib num2date and date2num
# These routines seem to give more consistent conversion.
# The original ones sometimes use a different base year, e.g. 1970. TODO: Figure out why and whether these routines are safe.
def num2date(mpltime):
    if np.ndarray in mpltime.__class__.__mro__:
        out = np.empty(len(mpltime), dtype='O')
        for idx, val in enumerate(mpltime.flat):
            out[idx] = num2date(val)
        out.shape = mpltime.shape
        return out
    return datetime.datetime.fromordinal(int(mpltime)) + datetime.timedelta(days=mpltime % 1)


def date2num(dt):
    if isinstance(dt, np.ndarray):
        if dt.dtype.name.startswith('datetime64'):
            dt = dt.astype('O')
        out = np.empty(len(dt), dtype=np.float64)
        for idx, val in enumerate(dt.flat):
            out[idx] = date2num(val)
        out.shape = dt.shape
        return out
    return (dt.toordinal() +
            (((dt.microsecond / 1e6 +
               dt.second) / 60 +
              dt.minute) / 60 +
             dt.hour) / 24)


def mat_able(d):
    """ Convert a data structure `d` to a format saveable by scipy.io.savemat """
    # print(type(d))

    maxlen = 63  # max length for Matlab field names: should be <31 characters (<63 for MATLAB 7.6+)

    def mat_field(f):
        # Modify a dict key to comply with Matlab variable naming rules
        mf = str(f).replace('.', '_').replace('-', '_') #TODO Any more substitutions needed?
        # should start with a letter
        if mf[0].isdigit():
            mf = 'NUMBER_' + mf
        # clip long names: cut out the middle and insert _CLIP_
        if len(mf) > maxlen:
            clipstr = '_CLIP_'
            ncut = (maxlen - len(clipstr)) / 2
            if clipstr in mf:
                # previously clipped string
                mfparts = mf.split(clipstr)
                mf = mfparts[0][:int(np.floor(ncut))] + clipstr + mfparts[1][-int(np.ceil(ncut)):]
            else:
                mf = mf[:int(np.floor(ncut))] + clipstr + mf[-int(np.ceil(ncut)):]
        return mf

    def undup(f):
        # Append numbers to duplicate strings in a list
        if len(f) != len(set(f)):
            counts = Counter(f)
            mf = [key if i == 0 else mat_field(key + '_' + str(i+1))
                  for key in set(f)
                  for i in range(counts[key])]
        else:
            mf = f
        return mf

    if isinstance(d, dict):
        # modify to make valid matlab field names
        mfields = [mat_field(key) for key in d.keys()]
        # check for duplicated field names
        mfields = undup(mfields)
        # if some names became too long, use numbered entries
        if np.any([len(k) for k in mfields]) > maxlen:
            mfields = ['field_{}'.format(k) for k in range(len(mfields))]
        dd = {key: value for key, value in zip(mfields, d.values())}
        # convert each item recursively
        for key, value in dd.items():
            # print('dict', key, type(value))  ###########
            dd[key] = mat_able(value)
        # add dict with original keys
        if set(mfields) != set(d.keys()):
            dd['ORIGINAL_KEYS'] = {key: value for key, value in zip(mfields, d.keys())}
    elif isinstance(d, datetime.datetime):
        dd = date2num(d) + 366
    elif isinstance(d, datetime.timedelta):
        dd = d.total_seconds()/3600/24  # to float days
    elif isinstance(d, np.ndarray):
        if d.dtype == 'O' and d.size > 0 and isinstance(d[0], datetime.datetime):
            dd = date2num(d) + 366
        else:
            dd = d
    elif isinstance(d, float):
        dd = d
    elif isinstance(d, np.float32):
        dd = d
    elif isinstance(d, np.int32):
        dd = d
    elif isinstance(d, np.int64):
        dd = d
    elif isinstance(d, int):
        dd = d
    elif isinstance(d, np.complex128):
        dd = d
    elif isinstance(d, complex):
        dd = d
    elif isinstance(d, str):
        dd = d
    elif isinstance(d, tuple):
        # print(type(d))###########
        dd = [mat_able(value) for value in d]
    elif isinstance(d, list):
        # print(type(d))  ###########
        dd = [mat_able(value) for value in d]
    else:
        # generic conversion: works for custom classes
        dd = {key: getattr(d, key) for key in dir(d) if
              not key.startswith("__") and not key.endswith("__")}
        # convert each item recursively
        for key, value in dd.items():
            # print('generic ', key, type(value))  ###########
            dd[key] = mat_able(value)

    return dd


def save_pickle_and_mat(opt, filepath, result):
    """ Generic routine to save analysis results """
    if opt['analysis']['save_pickle']:
        with open(filepath, 'wb') as fid:
            pickle.dump(result, fid, protocol=2)
    if opt['analysis']['save_mat']:
        # `long_field_names = True` assumes MATLAB 7.6+
        scipy.io.savemat(filepath.replace(".pickle", ".mat"),
                         {'class4': mat_able(result)},
                         long_field_names=True,
                         oned_as='column')

