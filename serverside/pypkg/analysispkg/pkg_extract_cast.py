import datetime
import functools
import multiprocessing as mp
import netCDF4 as nc
import numpy as np
import os
import time
import tqdm

from analysispkg import pkg_data
from analysispkg import pkg_eos
from analysispkg import pkg_extract
from analysispkg import pkg_grid
from analysispkg import pkg_interp
from analysispkg import pkg_obs
from analysispkg import pkg_utils


@pkg_utils.time_it
def extract_model(opt): # replaces model read in getMetrics
    # {
    #     "filename": "/home/sdfo600/DATA/OBS-WEST/CTD/20220210/101_CastCTD_2021-04-07_0056h.nc",
    #     "code": "101_CastCTD_2021-04-07_0056h",
    #     "lon": -124.75166320800781,
    #     "lat": 48.42333221435547,
    #     "time": [
    #         "2021-04-07T00:56:21",
    #         "2021-04-07T00:56:21"
    #     ]
    # },
    instr_type = "CTD"
    
    pattern   = opt['extract'][instr_type]['pattern']
    timename  = opt['extract'][instr_type]['time_var']

    mtimes,mfiles,mtimeind = pkg_data.nemo_files_and_times(opt, pattern, timename)

    # check if we have data to calculate stretched layer depths
    with nc.Dataset(mfiles[0]) as mCast:
        # layer depths from the stretched layer thicknesses, e3t
        if 'e3t' in mCast.variables:
            # e3t is recorded in ts files
            sshtimes, sshfiles, sshtimeind = None, None, None
        else:
            # no e3t in ts files; will have to get it from ssh files
            if 'ssh' in opt['extract'][instr_type].keys() and opt['extract'][instr_type]['ssh'] is not None:
                ssh_pattern = opt['extract'][instr_type]['ssh_pattern']
                sshtimes, sshfiles, sshtimeind = pkg_data.nemo_files_and_times(opt,ssh_pattern,timename)
                # check if ssh is recorded at same times as ts and chunked in same way
                print("len mtimes: ", len(mtimes))
                print("len sshtimes: ", len(sshtimes))              
                if len(mtimes) != len(sshtimes) or np.any(mtimes != sshtimes) or np.any(mtimeind != sshtimeind):
                    raise RuntimeError('SSH times must correspond to TS times in model output.')
                # TODO Also check if times are both instantaneous or both averaged?
            else:
                print('Neither `e3t` nor ssh files are available. No VVL stretching will be applied.')
                sshtimes, sshfiles, sshtimeind = 0, 0, 0  # use 0's as an indicator for no VVL stretching

    mtimes = pkg_data.apply_time_offset_seconds(mtimes, opt['extract'][instr_type]['time_offset'])

    # requested time frame for analysis
    mstart, mend = pkg_extract.get_period_bounds(opt)

    # Load observations index and filter by location, time, and exclude list
    observations = pkg_obs.observations_load_and_filter(opt, instr_type, mstart, mend)

    os.makedirs(opt['dir_extract'], exist_ok=True)
    logname = instr_type+'_extract_error'
    logname = logname + '_'+opt['runid'] + '.log'
    logfile = os.path.join(opt['dir_extract'],logname)

    if opt['parallel']:
        with mp.Pool(opt['nproc']) as pool:
            #results = [pool.apply_async(try_to_fetch_cast,
            #                            args=(opt, mfiles, mtimes, mtimeind, cast, sshfiles))
            #            for cast in observations]
                        
            # change made by MD August 2023
            results = []
            for cast in observations:
              i=np.searchsorted(mtimes, cast['time'][0])
              i1=max(0,i-3)
              i2=min(i+3,len(mfiles))
              results += [pool.apply_async(try_to_fetch_cast,
                                           args=(opt, 
                                                 mfiles[i1:i2], 
                                                 mtimes[i1:i2], 
                                                 mtimeind[i1:i2], 
                                                 cast, 
                                                 sshfiles[i1:i2]))
              ]
            
            
            errs = [result.get() for result in tqdm.tqdm(results)]
    else:
        errs = [try_to_fetch_cast(opt, mfiles, mtimes, mtimeind, cast, sshfiles) for cast in tqdm.tqdm(observations)]

    # write log file
    with open(logfile,'w') as f:
        f.write('Error messages from CTD cast analysis run at ' + str(datetime.datetime.utcnow()) + ':\n\n')
        nmissing = 0
        for e in errs:
            if e is not None:
                nmissing += 1
                f.write(e)
    print('Finished extraction. Number of casts failed to extract: {}'.format(nmissing))
    print('Log file: {}'.format(logfile))


@functools.lru_cache(maxsize=5)
def getgrid(file_coord, file_bathy, file_mesh):
    lonm, latm = pkg_grid.tgrid(file_coord)
    bathy = pkg_grid.bathymetry(file_bathy)
    tmaskm = pkg_grid.tmask(file_mesh)
    with nc.Dataset(file_mesh) as ncmesh:
        depm = ncmesh['gdept_0'][0,...] # 3D layer depths
        
    return lonm,latm,depm,tmaskm,bathy


def try_to_fetch_cast(opt,data_files,data_times,mtimeind,cast,sshfiles):
    """ Extracts a single cast """
    lonm, latm, depm, tmaskm, bathy = getgrid(opt['file_coord'],opt['file_bathy'],opt['file_mesh'])
    try:
        fetch_model_data(opt,data_files,data_times,mtimeind,cast,lonm,latm,depm,tmaskm,sshfiles,bathy)
        return None
    except Exception as e:
        errstr = 'Couldn\'t process cast ' + \
                    cast['filename'] + ' because:\n        ' + str(e) + '\n\n'
        return errstr


def fetch_model_data(opt,data_files,data_times,mtimeind,cast,lonm,latm,depm,tmaskm,sshfiles,bathy):
    """
    cast : obs dict from index
    period : model output interval, hrs
    """
    
    instr_type = "CTD"
    
    temp_name = opt['extract'][instr_type]['temperature']
    sal_name = opt['extract'][instr_type]['salinity']
    if 'ssh' in opt['extract'][instr_type].keys() and opt['extract'][instr_type]['ssh'] is not None:
        ssh_name = opt['extract'][instr_type]['ssh']
    else:
        ssh_name = 0 # use 0's as an indicator for no VVL stretching
    period = (data_times[1] - data_times[0]).total_seconds()/3600 # model write interval in hours

    #TODO Horizontal interpolation instead of nearest
    # Get interpolation weights
    # ws = pkg_interp.get_weights(lonm,latm,cast['lon'],cast['lat'],tmaskm)
    # if not ws['w'].any():
    #     # Fallback to nearest point if lininterp fails
    #     ws = pkg_interp.get_weights_nearest(cast['lon'],cast['lat'],opt['file_coord'],opt['file_bathy'])
    # i = slice(ws['i'][0],ws['i'][1]+1)
    # j = slice(ws['j'][0],ws['j'][1]+1)
    ws = pkg_interp.get_weights_nearest(cast['lon'], cast['lat'], opt['file_coord'], opt['file_bathy'],
                                        max_boxes=1, extrap=False)
    i = ws['i'][0]
    j = ws['j'][0]


    if type(tmaskm) is np.ma.masked_array:
        tmaskm = tmaskm.filled()
    if np.isnan(i) or tmaskm[0, j, i] == 0:
        # get_weights_nearest takes mask into account, so mask check is likely redundant here
        print('On land, skipping: ', cast['filename'])
        return

    # mLon, mLat, mDep = lonm[j,i], latm[j,i], depm[:,j,i]
    mLon, mLat = lonm[j, i], latm[j, i]
    water_depth = bathy[j, i]
    
    # #TODO Why do we trim to obs depths?
    # cD = np.where(np.logical_and(mDep >= oDep.min(),mDep <= oDep.max()))[0]
    # mDep = mDep[cD]

    # # find model time steps to interpolate to the obs time
    # from matplotlib.dates import date2num
    # i1, i2, w1, w2 = pkg_interp.interp_weights_1d(date2num(cast['time'][0]),
    #                                               date2num(data_times))

    # find file and time index for the nearest time
    tdif = np.abs(np.asarray(data_times) - cast['time'][0])
    ind = np.argmin(tdif)
    mFile = data_files[ind]
    mTime = data_times[ind]
    #return no data if the 
    if tdif[ind].total_seconds()/3600. > 2. * period:
        raise ValueError('No model data available within 2 output intervals of cast tiime.')
    
    # read model values and convert as needed
    with nc.Dataset(mFile) as mCast:
        if temp_name not in mCast.variables:
            raise NameError('No variable {} found in {}'.format(temp_name, mFile))
        if sal_name not in mCast.variables:
            raise NameError('No variable {} found in {}'.format(sal_name, mFile))
        flag_CT, flag_SR, flag_SA = pkg_eos.eoscheck(mCast, temp_name, sal_name)
        pTemp_m = mCast[temp_name][mtimeind[ind],:,j,i].filled()
        pSal_m  = mCast[sal_name][mtimeind[ind],:,j,i].filled()
        # layer depths from the stretched layer thicknesses, e3t

        if 'temperature_in_Kelvins' in opt['extract']['CTD'].keys() and opt['extract']['CTD']['temperature_in_Kelvins']:
            # print('Converting T from K to C')
            pTemp_m -= 273.15  # K to C
        # HACK: Remove 16th layer from files converted from rpn for CIOPS eval
        if 'layers_to_delete' in opt['extract']['CTD'].keys():
            # print('before remove: pTemp_m.shape',pTemp_m.shape)
            pTemp_m = np.delete(pTemp_m, opt['extract']['CTD']['layers_to_delete'])
            pSal_m = np.delete(pSal_m, opt['extract']['CTD']['layers_to_delete'])
            # print('after remove: pTemp_m.shape', pTemp_m.shape)

    # get layer depths
    if sshfiles == 0:
        
        # neither `e3t` nor ssh files are available, fallback to no VVL stretching
        mDep = depm[:, j, i]
        if depm.ndim == 1:
            mDep = depm
        else:
            mDep = depm[:, j, i]

    elif sshfiles is None:
        with nc.Dataset(mFile) as mCast:
            e3 = mCast['e3t'][mtimeind[ind], :, j, i].filled()
        mDep = np.concatenate(([0], np.cumsum(e3[:-1]))) + e3 * .5
    else:

        # NOTE: Assuming same times and chunking for ssh files as for ts files
        with nc.Dataset(sshfiles[ind]) as ncssh:
            if ncssh[ssh_name].ndim == 3:  # normal situation
                ssh_m = ncssh[ssh_name][mtimeind[ind], j, i].filled()
            elif ncssh[ssh_name].ndim == 4:  # files converted from rpn (ciops eval)
                ssh_m = ncssh[ssh_name][mtimeind[ind], 0, j, i].filled()
            else:
                raise ValueError('SSH variable should have 3 or 4 dimensions, but has shape {}'.format(ncssh[ssh_name]))
        mDep = depm[:,j,i] * (1 + ssh_m/water_depth)  # vvl stretching

    # apply mask
    # with nc.Dataset(opt['file_mesh']) as ncm:
    #     maskji = ncm['tmask'][0,:,j,i].filled()
    maskji = tmaskm[:, j, i]
    ival = np.where(maskji==1)[0]
    pTemp_m = pTemp_m[ival]
    pSal_m = pSal_m[ival]
    mDep = mDep[ival]

    # keep only t,s where both of them have valid values  # MD: Is this neccessary? Does NEMO model ever produce s<0 or t<-2?
    cT = np.where(np.logical_and(np.isfinite(pTemp_m),pTemp_m > -2.))[0]
    cS = np.where(np.logical_and(np.isfinite(pSal_m),pSal_m > 0.))[0]
    c = np.intersect1d(cT,cS)
    pTemp_m = pTemp_m[c]
    pSal_m = pSal_m[c]
    mDep = mDep[c]

    # Handle TEOS-10 to EOS-80
    if flag_CT and flag_SR:
        pTemp_m, pSal_m = pkg_eos.CTSR_to_PTSP(pTemp_m, pSal_m, mDep, mLon, mLat)
    if flag_CT and flag_SA:
        pTemp_m, pSal_m = pkg_eos.CTSA_to_PTSP(pTemp_m, pSal_m, mDep, mLon, mLat)

    mod = {'z':mDep, 'pTemp':pTemp_m, 'salinity':pSal_m,
         'water_depth':water_depth, 'time':mTime, 'lat':mLat,'lon':mLon}
    #TODO Extract water depth from bathy file?

    # Save model extraction
    pfile = pkg_data.save_extraction(opt, cast, mod)


    # Confirm that we can reload the pickle
    try:
        mod, pfile = pkg_data.load_extraction(opt, cast)
    except Exception as err:
        print("Failed to reload pickle {}\n{}".format(pfile,err))

