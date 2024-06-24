#  Python utility to generate a YAML index of each instrument in an observations archive

import argparse
import datetime
import glob
import json
import multiprocessing as mp
import netCDF4 as nc
import numpy as np
import os
import time

from analysispkg import pkg_obs


def read_time_frame(ncid, f,  timename='time', calendar='standard'):
    """
    Utility: reads the specified time variable from a netcdf file, return the min,max values as datetime objects.

    It is expensive to run nc.num2date on long time series, and for indexing we only need the min and max dates. So
    here we only convert those two numbers.
    """
    if timename not in ncid.variables:
      timename = 't'
      if 't' not in ncid.variables: # GO added
        raise RuntimeError('Time variable {} not found in file {}.'.format(timename,f))

    time_raw = ncid.variables[timename][:]
    units = ncid.variables[timename].units

    time_frame_raw = [np.min(time_raw).astype(int), np.max(time_raw).astype(int)]

    # converts the numbers from the file to datetime objects
    time_frame = nc.num2date(time_frame_raw, units, calendar=calendar)

    # convert to datetime object
    time_frame = [datetime.datetime.fromtimestamp(t._to_real_datetime().timestamp()) for t in time_frame]

    return time_frame


def index_one_file(path,instrument,f):
    # Open the file, get lon and lat
    with nc.Dataset(f) as ncf:
        try:
            obs_lon = ncf.variables["longitude"][:]
            obs_lat = ncf.variables["latitude"][:]
            time_frame = read_time_frame(ncf,f)
            id_code = pkg_obs.get_station_label(ncf)
        except Exception as err:
            print("Failed to scan file {}".format(f))
            print(err)
            raise
            #continue

    lon = float(obs_lon)
    lat = float(obs_lat)
    meta = {
        "filename": os.path.relpath(f, start=path),
        "code": id_code,
        "lon": lon,
        "lat": lat,
        "time": time_frame
    }
    return meta


def generate_index_list(path, instr, mpool):
    instdir = os.path.join(path, instr)

    if not os.path.exists(instdir):
        # print("Did not find a directory for {}, skipping".format(instr))
        return None

    print("Indexing {}".format(instr))
    search_string = os.path.join(instdir, "**/", "*.nc")
    files = glob.glob(search_string, recursive=True)
    files.sort()

    # Index each file
    if mpool is not None:
        iterator = [mpool.apply_async(index_one_file, args=(path, instr, f,)) for f in files]
        index_list = [r.get() for r in iterator]
    else:
        index_list = [index_one_file(path, instr, f) for f in files]
    return index_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Index observations archive")
    parser.add_argument("path", type=str,
                        help="path to scan observations archive.")
    parser.add_argument("--indexpath", type=str,
                        help="Path to write index files, if not provided we write to INDEX under the observations archive path.")
    parser.add_argument("-j", type=int,
                        help="Number of processors to use for parallelized indexing.")
    args = parser.parse_args()

    # Types of instruments that we can index
    types = ['CTD', 'LH', 'SST']

    # Decide default pool size based on hostname
    if "inter-dfo" in os.uname().nodename:
        nproc = 40  # gpsc7 interactive node
    elif "ib14be-" in os.uname().nodename or "vis-" in os.uname().nodename:
        nproc = 64  # gpsc7 compute node
    else:
        nproc = 1

    # Override pool size if requested on cmdline
    if args.j is not None:
        if args.j > 0:
            nproc = args.j

    # By default we write indices to INDEX subdir, unless path provided at cmdline
    if args.indexpath is None:
        indexpath=os.path.join(args.path,'INDEX')
    else:
        indexpath=args.indexpath
    os.makedirs(indexpath,exist_ok=True)

    t0=time.time()

    if nproc > 1:
        pool = mp.Pool(processes=nproc)
    else:
        pool = None

    # Loop over instruments and do the indexing
    for instrument in types:
        location_list = generate_index_list(args.path, instrument, pool)
        if location_list is None:
            print("No data for {} in {}, skipping".format(instrument, args.path))
            continue

        # Store result in a JSON file
        indexfile = os.path.join(indexpath, instrument + ".json")
        with open(indexfile, "w") as outfile:
            json.dump(location_list,outfile, default=datetime.datetime.isoformat, indent=2)

        print("Indexing {} complete: {} observations indexed".format(instrument,len(location_list)))

    if nproc > 1:
        pool.close()
        pool.join()

    print("obsindex.py completed indexing in {:.3} seconds\n\n".format(time.time()-t0))
