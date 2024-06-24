import datetime
import multiprocessing as mp
import netCDF4 as nc
import numpy as np
import tqdm
from scipy.spatial import cKDTree

from analysispkg import pkg_data
from analysispkg import pkg_eos
from analysispkg import pkg_extract
from analysispkg import pkg_grid
from analysispkg import pkg_interp
from analysispkg import pkg_obs
from analysispkg import pkg_signal


# made an edit for SST and SSS for LH and SST to be grabbed from z = 1 instead of 0
go_changes = 'off' 

def extract_model_timeseries(opt,instr_type):
    """ Extract model data corresponding to each observation

    - Time is not interpolated (we keep the model times)
    - If TG or SST, we attempt to horizontally interpolate first, and fall back to nearest grid point
    - Extracted result is saved as a dictionary in a pickle file

    Parameters
    ----------
    opt : dict
        Options dict.
    instr_type : 'SST' | 'LH'
        Type of instrument.
    """
    print('Begin {} time series extraction for {}-{}'.format(instr_type, opt['config'],opt['runid']))

    glamt, gphit = pkg_grid.tgrid(opt['file_coord'])
    tmaskm = pkg_grid.tmask(opt['file_mesh'])
    zm = pkg_grid.nav_lev(opt['file_mesh'])

    # requested time frame for analysis
    mstart, mend = pkg_extract.get_period_bounds(opt)

    print('Loading obs database...')
    # Load observations index and filter by location, time, and exclude list
    observations = pkg_obs.observations_load_and_filter(opt, instr_type, mstart, mend)

    if len(observations) == 0:
        print("No observations; exiting")
        return
    # Load obs data and truncate to model time bounds
    if opt['parallel']:
        with mp.Pool(processes=min(opt['nproc'], len(observations))) as pool:
            iterator = [pool.apply_async(pkg_obs.load_observation,
                                         args=(observation['filename'], instr_type, opt))
                        for observation in observations]
            obsdata = [r.get() for r in iterator]
    else:
        obsdata = [pkg_obs.load_observation(observation['filename'], instr_type, opt)
                   for observation in observations]

    obsdata = [pkg_obs.truncate_obs_timeseries(obs, mstart, mend) for obs in obsdata]
    obsdata = list(filter(None, obsdata))  # drop any None entries


    for obs in obsdata:
        print('Calculating weights for station ', obs['station'])
        # Get interpolation weights
        ws = pkg_interp.get_weights(glamt,gphit,obs['lon'],obs['lat'],tmaskm)
        if not ws['w'].any():
            # Fallback to nearest point if lininterp fails
            ws = pkg_interp.get_weights_nearest(obs['lon'],obs['lat'],opt['file_coord'],opt['file_bathy'])

        # Add weights to obs
        obs['ws'] = ws
        if opt['verbose'] >= 2:
            print("extract_one_observation: stn {} python indices {},{}".format(obs['station'], ws['i'], ws['j']))

    obsdata = [obs for obs in obsdata if 'ws' in obs.keys()]
    print("{} observations to extract after weights calculation".format(len(obsdata)))

    if len(obsdata) == 0:
        print('No observations extractable')
        return

    # Extraction time window
    tstart, tend = mstart, mend
    tmin = min([min(obs['time']) for obs in obsdata])
    tmax = max([max(obs['time']) for obs in obsdata])
    tstart, tend = max(tstart,tmin), min(tend, tmax)
    print("Extraction time window, start:", tstart, "end:", tend)

    # Subset the obsdata list to minimize size of obsdata passed to parallel pool
    keys = ['station', 'ws', 'lon', 'lat', 'inst_depth', 'shortname']
    obslist = [{key: value for key, value in obs.items() if key in keys} for obs in obsdata]

    # Extract
    print('Finding model files ...')
    pattern = opt['extract'][instr_type]['pattern']
    timename = opt['extract'][instr_type]['time_var']
    mfiles = pkg_data.list_source_files_sorted_in_range(opt, pattern, timename, tstart, tend)
    if len(mfiles) == 0:
        print("No model files found for pattern {}, can not proceed".format(pattern))
        return

    print('Begin extraction...')
    moddata = extract_observations(opt, instr_type, obslist, mfiles, tstart, tend)
    print("Extracted {} stations".format(len(moddata)))

    # Save model extractions
    for obs,mod in zip(obsdata,moddata):
        pfile = pkg_data.save_extraction(opt,obs,mod)
        print("Saved model extraction to {}".format(pfile))



def extract_observations(opt,instr_type,obslist,mfiles,tstart,tend, using_vvl=False, mfiles_vvl=None, bathy=None, zm=None ):
    """ Extract model time series corresponding to each instrument record.
    """

    time_offset = opt['extract'][instr_type]['time_offset']
    time_var = opt['extract'][instr_type]['time_var']

    if instr_type == 'SST':
        loader = load_data_sst
        tname = opt['extract'][instr_type]['temperature']
        argslist = [(mfile, tstart, tend, time_offset, tname, obslist, time_var) for mfile in mfiles]

    if instr_type == 'LH':
        loader = load_data_lh
        tname = opt['extract'][instr_type]['temperature']
        sname = opt['extract'][instr_type]['salinity']
        argslist = [(mfile, tstart, tend, time_offset, tname, sname, obslist, time_var) for mfile in mfiles]

    if opt['parallel']:
        # Read mfiles in parallel
        nproc=min(opt['nproc'], len(mfiles))
        with mp.Pool(processes=nproc) as pool:
            iterator = [pool.apply_async(loader,args=args) for args in argslist]
            results = [r.get() for r in tqdm.tqdm(iterator)]

    else:
        results = []
        for args in tqdm.tqdm(argslist):
            results += [loader(*args)]

    # Filter out None entries; results is indexed as [mfile][time, [var][obs]]
    results = list(filter(lambda x: x is not None, results))

    # Unpack so time is indexed by [mfile] and mdatas by [var][mfile][obs]
    times, *mdatas = zip(*results)

    # Concatenate the times into final time vector, ensure it's sorted ascending
    time = np.concatenate(times)
    isort = time.argsort()
    time = time[isort]


    mdata=[]
    for var in mdatas:
        # Transpose from [mfile][obs] to [obs][mfile]
        tmp = list(map(list, zip(*var)))

        # Concatenate the data in the [mfile] dimension
        tmp2 = [np.concatenate(x) for x in tmp]
        # Apply interpolation weights to reduce the extractions into single time series
        tmp3 = [pkg_interp.interp(x, obslist[i]['ws']['w']) for i, x in enumerate(tmp2)]
        # Apply sorting
        tmp4 = [x[isort] for x in tmp3]
        mdata += [tmp4]
    mdata = tuple(mdata)

    if instr_type == 'SST':
        sst, = mdata
        d = [{'temperature':sst[i], 'time':time} for i in range(len(sst))]

    if instr_type == 'LH':
        sst, sss = mdata
        d = [{'temperature':sst[i], 'salinity':sss[i], 'time':time} for i in range(len(sst))]

    for i in range(len(obslist)):
        d[i]['station'] = obslist[i]['station']
    return d


def load_data_sst(mf, tstart, tend, time_offset, tname, obslist, time_var):
    """ Read SST data for specified indices and time range from one model file.
    """
    with nc.Dataset(mf) as ncf:
        c,timem = pkg_extract.load_time(ncf, tstart, tend, time_offset=time_offset, time_var=time_var)
        if c is None:
            return None
        # Load entire variable
        if ncf[tname].ndim == 3:    # Loading from 2D SST output
            sst = ncf[tname][...]
        elif ncf[tname].ndim == 4:  # Loading from top layer of 3D SST output
            if go_changes == 'on':
              #print('changes on - drawing from z = 1 instead of 0 -GO')
              sst = ncf[tname][:, 1, ...]
            else:
              sst = ncf[tname][:, 0, ...]
    # Extract each observation
    sst_out = []
    for obs in obslist:
        i,j,_ = toslice(obs['ws'])
        sst_out += [sst[c,j,i]]
    return timem, sst_out


def load_data_lh(mf, tstart, tend, time_offset, tname, sname, obslist, time_var):
    """ Lighthouse data loader """
    t = load_data_sst(mf, tstart, tend, time_offset, tname, obslist, time_var)
    s = load_data_sst(mf, tstart, tend, time_offset, sname, obslist, time_var)
    if t is None or s is None:
        return None
    else:
        timemt, sst_out = t
        timems, sss_out = s
        return timemt, sst_out, sss_out


def get_bathy (obslist, bathy):
    bathys_out = []
    for obs in obslist:
        i,j,_ = toslice(obs['ws'])
        bathys_out += [pkg_interp.interp(bathy[j,i], obs['ws']['w'][0]) ]
    return bathys_out


def toslice(ws):
    i = slice(ws['i'][0],ws['i'][1]+1)
    j = slice(ws['j'][0],ws['j'][1]+1)
    k = slice(ws['k'][0],ws['k'][1]+1)
    return i,j,k
