import datetime
import functools
import glob
import json
import multiprocessing as mp
import netCDF4 as nc
import numpy as np
import os
import time
import tqdm

from analysispkg import pkg_data

"""
    scan() opens every netcdf file found under opt['src_output_dir'], loads the time vector info,
    and caches the result in a JSON file.

    To load the index:
    timedata = pkg_scan.load_netcdf_index(opt)

    To get the time vector:
    mtime = pkg_scan.get_time(timedata, filename, timevectorname)
"""

def scan(opt,force=False):
    """
    Load the time vectors from every netcdf and cache it as a dict in a JSON file:
     - If the cache does not exist, or if force=True, we scan all netcdf files and write a new index.
     - If the cache exists, we scan any unscanned files and write an updated index.
     - If cache exists and there are no unscanned files we do nothing.
    """
    root = opt['src_output_dir']
    allfiles = pkg_data.list_source_files(root, "*.nc")

    indexfile = indexfilename(opt)
    if not os.path.exists(indexfile) or force:
        timedata = {}
        files = allfiles
    elif os.path.exists(indexfile):
        timedata = load_netcdf_index(opt)
        scanned_files = [os.path.join(root, f) for f in timedata.keys()]
        files = list(set(allfiles) - set(scanned_files))
        print("Found {} netcdf files; {} already scanned, {} to scan".format(len(allfiles), len(scanned_files), len(files)))
    else:
        raise RuntimeError("Impossible situation in pkg_scan.scan()")

    nf = len(files)
    if nf == 0:
        print("No files to scan.")
        return

    print("Begin scanning {} netcdf files ...".format(len(files)))
    t0 = time.time()
    if opt['parallel']:
        with mp.Pool(processes=opt['nproc']) as pool:
            iterator = [pool.apply_async(get_timeinfo, args=(f,root)) for f in files]
            for x in tqdm.tqdm(iterator):
                file, timeinfo = x.get()
                timedata[file] = timeinfo
    else:
        for f in files:
            file, timeinfo = get_timeinfo(f,root)
            timedata[file] = timeinfo
    t1 = time.time()
    print('Scanned {} files in {:.3f}s '.format(nf,t1-t0))

    # Drop any files with zero time records
    dropfiles=[]
    for f in timedata.keys():
        for v in timedata[f].keys():
            if len(timedata[f][v]['data']) == 0:
                print("File {} var {} has length 0".format(f,v))
                dropfiles += [f]
    for f in set(dropfiles):
        print('Excluding {} from the index'.format(f))
        del timedata[f]

    t0 = time.time()
    save_netcdf_index(opt, timedata)
    t1 = time.time()
    print('Saved JSON index to {} in {:.3f}s '.format(indexfile,t1-t0))

    t0 = time.time()
    _ = load_netcdf_index(opt)
    t1 = time.time()
    print('Test reloading JSON index in {:.3f}s '.format(t1-t0))


def get_timeinfo(file,root):
    """ Loads the time vector and metadata """
    timeinfo = {}
    with nc.Dataset(file) as ncf:
        for var in set(ncf.variables.keys()).intersection(set(['time_counter','time_instant'])):
            timeinfo[var] = {
                'data': list(ncf[var][:].filled()),
                'meta': ncf[var].__dict__,
            }
    relfile = os.path.relpath(file, start=root)
    return relfile,timeinfo


def indexfilename(opt):
    return os.path.join(opt['dir_extract'], "NETCDF_INDEX.json")


def save_netcdf_index(opt,timedata):
    """ Save the netcdf index """
    indexfile = indexfilename(opt)
    indexpath, _ = os.path.split(indexfile)
    os.makedirs(indexpath, exist_ok=True)
    with open(indexfile, "w") as outfile:
        json.dump(timedata, outfile, indent=1)


@functools.lru_cache(maxsize=5)
def load_json_from_file(infile,root):
    with open(infile, "r") as fid:
        data = json.load(fid)
    # relpath to fullpath
    data2 = {os.path.join(root,k): v for k, v in data.items()}
    return data2


def load_netcdf_index(opt):
    """ Load the netcdf index
     """
    indexfile = indexfilename(opt)
    if not os.path.exists(indexfile):
        raise FileNotFoundError("Model netcdf index {} not found, please run scan.py!".format(indexfile))
    root = opt['src_output_dir']
    timedata = load_json_from_file(indexfile,root)
    return timedata


def get_time(timedata,file,timename):
    """
    Consult the timedata index and return the datetimes for a given file,timename
    """
    time_raw = timedata[file][timename]['data']
    units = timedata[file][timename]['meta']['units']
    calendar = timedata[file][timename]['meta'].get('calendar','standard')

    if 'seconds since' in units:
        # Convert reference time to netcdf4's datetime, then to datetime.datetime
        time00 = nc.num2date(0, units, calendar)
        time0 = datetime.datetime(time00.year, time00.month, time00.day, time00.hour, time00.minute, time00.second)
        # Create ufunc for converting seconds to dt.timedelta
        to_timedelta = np.frompyfunc(lambda s: datetime.timedelta(seconds=s),1,1)
        # Apply ufunc to time_raw and add reference time to get array of dt.datetime
        time_var = to_timedelta(time_raw) + time0
    else:  # old method
        time_var = nc.num2date(time_raw, units, calendar)
        # Convert to array of datetime.datetime
        time_var = np.array([datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in time_var])

    return time_var
