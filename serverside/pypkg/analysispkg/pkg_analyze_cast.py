import datetime as datetime
import multiprocessing as mp
import numpy as np
import os
import warnings

from analysispkg import pkg_data
from analysispkg import pkg_obs
from analysispkg import pkg_statistics
from analysispkg import pkg_utils

# required to get domain - GO 20230612
from matplotlib.path import Path as Path
# required to get tmask, e3t0 mask etc - GO 20230919
#import netCDF4 as nc
go_changes = "on" # switch for my changes - GO 20230612

def analyze_casts(opt):
    paths = {}
    paths['model_data_directories'] = [opt['src_output_dir']]
    paths['analysis_base'] = os.path.join(opt['dir_process'], 'CTD')
    paths['bounding_boxes'] = os.path.join(paths['analysis_base'],'bounding_Boxes_CTDAnalysis.pickle')
    paths['figure_output_base_folder'] = opt['dir_plots'] + '/CTD/'

    casename = opt['casename']
    runid = opt['runid']

    os.makedirs(paths['analysis_base'], exist_ok=True)

    if go_changes == 'on':
        depth_levels,sdomains = read_CTD_domain(opt['analysis']['CTD']['domain_file'])
    else:
        depth_levels,_ = read_CTD_domain(opt['analysis']['CTD']['domain_file'])
#        file_mesh = opt['file_mesh']
#        
#        with nc.Dataset(file_mesh) as mesh:
#      #     print(mesh.variables)
#          tmask=mesh.variables['tmask'][:] # 0's where depth lev exceeds
#          e3t0=mesh.variables['e3t_0'][:] # 'widths' or 'weights' for each depth lev
#          gdept_0=mesh.variables['gdept_0'][:] # depth levels (t; centre)
        

    #run validation metrics (this is only done once, then saved as a pickle file for further analysis and plot refinement)
    print('Calculating validation metrics against provided archive of CTD casts at ' +
          opt['obs']['root'] + ' for the run ' + casename + " ...")

    for period, (ana_start, ana_end) in opt['analysis']['periods'].items():
        print("Working on period {}".format(period))

        logfile = os.path.join(paths['analysis_base'], 'CTDcast_analysis_errorLog_{}.txt'.format(period))
        pfile = os.path.join(paths['analysis_base'],'CTDcast_metrics_{}.pickle'.format(period))

        if go_changes == 'on':
            metrics = getMetrics_GO(opt,depth_levels,sdomains,ana_start, ana_end, logfile)
        else:
            metrics = getMetrics(opt,depth_levels, ana_start, ana_end, logfile)



        print('Saving CTD cast skill metrics to ' + pfile)
        pkg_data.save_pickle_and_mat(opt, pfile, metrics)

        print('Saved CTD cast skill metrics to ' + pfile)
        print(' ')

    print('Finished.')


def read_CTD_domain(domain_file):
    try:
        data = pkg_utils.load_yaml(domain_file)
    except FileNotFoundError:
        print("WARNING:\n CTD domain_file {} not found as a relative or absolute path".format(domain_file))
        print(" We will try to find the domain file directly from config; this fallback will be removed in future cleanups")
        data = pkg_utils.load_config_yaml(domain_file)

    dep = np.asarray(data['depth']).astype(float)
    depth_levels = [[dep[i],dep[i+1]] for i in range(len(dep)-1)]

    coords = {}
    for c in data['polygon_coords'].keys():
        coords[c] = np.asarray(data['polygon_coords'][c])
        # make sure the polygon is closed
        if not np.all(coords[c][-1] == coords[c][0]):
            coords[c] = np.vstack([coords[c], coords[c][0, :]])

    return depth_levels,coords


def getSkill(obsData,modData,depth_levels):

    # Unpack data_classes.ctd classes to local variables
    pTemp_o,pSal_o,oDep = obsData['pTemp'], obsData['salinity'], obsData['z']
    time = obsData['time']
    pos = obsData['lat'], obsData['lon']

    pTemp_m,pSal_m,mDep = modData['pTemp'], modData['salinity'], modData['z']
    mTime = modData['time']
    mCoords = np.hstack([modData['lat'], modData['lon']])

    cDep = np.logical_or(mDep < oDep.min(),mDep > oDep.max())

    #interpolate observed data to model depth levels to create class4 data
    T_interp = np.interp(mDep,oDep,pTemp_o)
    S_interp = np.interp(mDep,oDep,pSal_o)
    T_interp[cDep] = np.nan
    S_interp[cDep] = np.nan
    if np.sum(np.isfinite(T_interp)) <= 2:
        return None, None, 0

    class4 = {'obs':{'T':T_interp,'S':S_interp},'model':{'T':pTemp_m,'S':pSal_m},'dep':mDep}

    labels = ['Full']
    levels = [[0.,mDep[-1]]]
    
    for lev in depth_levels:
        labels.append('%1i'%lev[0] + '->%1i'%lev[1])
        levels.append(lev)
    labels.append('>%1i'%depth_levels[-1][1])
    levels.append([depth_levels[-1][1],mDep[-1]])

    scores = {'bias':{},'rmse':{},'crmse':{},'skill1981':{},'ModelTime':mTime,'ModelLocation':mCoords,'ObsTime':time,'ObsLocation':pos}
    for lab,rang in zip(labels,levels):
        for metric in ['bias','rmse','crmse','skill1981']:
            scores[metric][lab] = {}
        for var in ['T','S']:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                bias, rmse, crmse, skill1981 = calcScores(var,rang,mDep,class4)               
                                
            scores['bias'][lab][var] = bias
            scores['rmse'][lab][var] = rmse
            scores['crmse'][lab][var] = crmse
            scores['skill1981'][lab][var] = skill1981

    return scores, class4, 1

# duplicate of above but includes mean and std dev of mod / obs and subdomain - GO 20230612
def getSkill_GO(obsData,modData,depth_levels,sdomains):
  
    # Unpack data_classes.ctd classes to local variables
    pTemp_o,pSal_o,oDep = obsData['pTemp'], obsData['salinity'], obsData['z']
    time = obsData['time']
    pos = obsData['lat'], obsData['lon']
    
    pTemp_m,pSal_m,mDep = modData['pTemp'], modData['salinity'], modData['z']
    mTime = modData['time']
    mCoords = np.hstack([modData['lat'], modData['lon']])
    
    cDep = np.logical_or(mDep < oDep.min(),mDep > oDep.max())
    
    # get subdomain - GO 20230612
    sdcast = ""
    for sd in sdomains.keys():
      if Path(sdomains[sd]).contains_point(pos):
        sdcast = sd

    #interpolate observed data to model depth levels to create class4 data
    T_interp = np.interp(mDep,oDep,pTemp_o)
    S_interp = np.interp(mDep,oDep,pSal_o)
    T_interp[cDep] = np.nan
    S_interp[cDep] = np.nan
    if np.sum(np.isfinite(T_interp)) <= 2:
        return None, None, 0

    class4 = {'obs':{'T':T_interp,'S':S_interp},'model':{'T':pTemp_m,'S':pSal_m},'dep':mDep}

    labels = ['Full']
    levels = [[0.,mDep[-1]]]
    
    for lev in depth_levels:
        labels.append('%1i'%lev[0] + '->%1i'%lev[1])
        levels.append(lev)
    labels.append('>%1i'%depth_levels[-1][1])
    levels.append([depth_levels[-1][1],mDep[-1]])

    scores = {'bias':{},'rmse':{},'crmse':{},'skill1981':{},'mean_obs':{},'mean_mod':{},'std_obs':{},'std_mod':{},'ModelTime':mTime,'ModelLocation':mCoords,'ObsTime':time,'ObsLocation':pos, 'SubDomain':sdcast}
    
    for lab,rang in zip(labels,levels):
        for metric in ['bias','rmse','crmse','skill1981','mean_obs','mean_mod','std_obs','std_mod']:
            scores[metric][lab] = {}
        for var in ['T','S']:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                
                bias, rmse, crmse, skill1981, mean_obs, mean_mod, std_obs, std_mod = calcScores_GO(var,rang,mDep,class4)               
                                
            scores['bias'][lab][var] = bias
            scores['rmse'][lab][var] = rmse
            scores['crmse'][lab][var] = crmse
            scores['skill1981'][lab][var] = skill1981
            scores['mean_obs'][lab][var] = mean_obs
            scores['mean_mod'][lab][var] = mean_mod
            scores['std_obs'][lab][var] = std_obs
            scores['std_mod'][lab][var] = std_mod

    return scores, class4, 1

def calcScores(var,dRange,mDep,class4):
    #return all NaN's if specified depth range is deeper than model bathymetry at cast location
    if dRange[0] > mDep[-1]:
        return np.nan, np.nan, np.nan, np.nan

    c = np.where(np.logical_and(mDep >= dRange[0],mDep <= dRange[1]))[0]

    mod = class4['model'][var][c]
    obs = class4['obs'][var][c]

    bias = pkg_statistics.bias(obs,mod)
    rmse = pkg_statistics.rmse(obs,mod)
    crmse = pkg_statistics.crmse(obs, mod)
    skill1981 = pkg_statistics.Willmott1981(obs,mod)

    return bias,rmse,crmse,skill1981

# duplicate of calcScores above but w/ simple mean and std included - GO 20230612
# GO 20230919 - note the simple mean incorrect. should account for class4 on model stretched Z
def calcScores_GO(var, dRange, mDep, class4):
    #return all NaN's if specified depth range is deeper than model bathymetry at cast location
    if dRange[0] > mDep[-1]:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    c = np.where(np.logical_and(mDep >= dRange[0],mDep <= dRange[1]))[0]

    mod = class4['model'][var][c]
    obs = class4['obs'][var][c]

    if (mod.shape == (1,)) & (obs.shape == (1,)):
      return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      
    bias = pkg_statistics.bias(obs,mod)
    rmse = pkg_statistics.rmse(obs,mod)
    crmse = pkg_statistics.crmse(obs, mod)
    skill1981 = pkg_statistics.Willmott1981(obs,mod)
    
    # GO 20230612 - fixes to calc deptht weighted vert avg and stdv
    # TODO: using mDep to infer e3t - actually get e3t from mesh in future
    e3t = np.zeros_like(mDep)
    e3t[0] = mDep[0] * 2 # shallowest depth bin width
    for i in range(1, len(mDep)):
      e3t[i] = (mDep[i] - (mDep[i - 1] + (e3t[i - 1] / 2))) * 2
    # truncate e3t to match data shape
    indices = np.where((mDep >= dRange[0]) & (mDep <= dRange[1]))[0]
    e3t = e3t[indices[0]:indices[-1]+1]
    mDep2 = mDep[indices[0]:indices[-1]+1]
    
#    mean_obs = pkg_statistics.nanmean_simple(obs)
#    mean_mod = pkg_statistics.nanmean_simple(mod)

    # need to have a 'nan mask' based on obs nans and pass that
    # obsnanmask
    nan_indices = np.isnan(obs)
    
    mod[nan_indices] = np.nan
    e3t[nan_indices] = np.nan
#    print('mod after nan ', mod)
    #print('e3t after nan ', e3t)
    
    mean_obs = pkg_statistics.nanmean_vvl(obs, mDep2, e3t) # new functions added and tested
    mean_mod = pkg_statistics.nanmean_vvl(mod, mDep2, e3t)
    
#    print('the obs ', obs)
#    print('the c ', c)
#    print('mean obs old way ', mean_obs)
#    print('mean mod new way ', mean_mod)
#    print('mean obs new way ', mean_obs)
#    print('mean mod new way ', mean_mod)
    
    #std_obs = pkg_statistics.nanstd_simple(obs)
    #std_mod = pkg_statistics.nanstd_simple(mod)
    std_obs = pkg_statistics.nanstd_vvl(obs, mDep2, e3t)
    std_mod = pkg_statistics.nanstd_vvl(mod, mDep2, e3t)
    
    return bias,rmse,crmse,skill1981,mean_obs,mean_mod,std_obs,std_mod

#function to calculate performance metrics for 9 grid points in immediate vicinity of cast location
#currently unused and not thoroughly tested, but could be implemented for future functionality
def calcScores_multiPoint(var,dRange,mDep,class4):
    if dRange[0] > mDep[-1]:
        return np.nan, np.nan, np.nan

    c = np.where(np.logical_and(mDep >= dRange[0],mDep <= dRange[1]))[0]

    skill = np.zeros([3,3]); bias = np.zeros([3,3]); rmse = np.zeros([3,3])
    for i in range(class4['model'][var].shape[1]):
        for j in range(class4['model'][var].shape[2]):
            n = np.sum((class4['model'][var][c,i,j] - class4['obs'][var][c])**2.)
            dM = np.abs(class4['model'][var][c,i,j] - np.nanmean(class4['model'][var][c,i,j]))
            dO = np.abs(class4['obs'][var][c] - np.nanmean(class4['obs'][var][c]))
            d = np.sum((dM + dO)**2.)

            skill[i,j] = np.max([0.,1. - n/d])
            bias[i,j] = np.nanmean(class4['model'][var][c,i,j] - class4['obs'][var][c])
            rmse[i,j] = np.sqrt(np.nanmean((class4['model'][var][c,i,j] - class4['obs'][var][c])**2.))

    select = skill / (bias*rmse)
    if (select == 0.).all():
        select = np.ones_like(skill) / (bias*rmse)

    c = np.where(select.flatten() == np.nanmax(select))[0]
    if len(c) == 0:
        c = [5]
    sBias = bias.flatten()[c]
    sRMSE = rmse.flatten()[c]
    sSkill = skill.flatten()[c]

    return float(sBias[0]),float(sRMSE[0]),float(sSkill[0])

def getMetrics(opt,depth_levels,ana_start,ana_end,logfile):
    casename = opt['casename']
    metrics = {casename: {'scores':{},'class4':{}}}
    missedFiles = 0

    # Load observations index and filter by location, time, and exclude list
    listCasts = pkg_obs.observations_load_and_filter(opt, "CTD", ana_start, ana_end)
    totalFiles = len(listCasts)
    if totalFiles == 0:
        print('No casts to process')
        return

    with open(logfile,'w') as f:
        f.write('Error messages from CTD cast analysis run at ' + str(datetime.datetime.utcnow()) + ':\n\n')

        print("Begin calculating metrics for {} casts ...".format(totalFiles))
        if opt['parallel']:
            with mp.Pool(processes=min(opt['nproc'], len(listCasts))) as pool:
                iterator = [pool.apply_async(process_a_cast, args=(opt, cast, depth_levels)) for cast in listCasts]
                results = [r.get() for r in iterator]
        else:
            results = [process_a_cast(opt, cast, depth_levels) for cast in listCasts]

        # Put the results in the output structure
        for result in results:
            scores, class4s, name, e = result
            if e is not None:
                if e != 'single pt':
                    missedFiles += 1
                    f.write(str(missedFiles) + '. Couldn\'t process cast ' + name + ' for ' + casename + ' because:\n        ' + str(e) + '\n\n')
            else:
                metrics[casename]['scores'][name],metrics[casename]['class4'][name] = scores, class4s

    print('Finished processing {} casts.'.format(totalFiles))
    if missedFiles > 0:
        print('{} casts could not be analyzed, see error log: {}'.format(missedFiles,logfile))

    return metrics

# only change to above getMetrics is sdomains is passed and process_a_cast_GO is called instead
def getMetrics_GO(opt,depth_levels,sdomains,ana_start,ana_end,logfile):
     
    print("in getMetrics_GO")

    casename = opt['casename']
    metrics = {casename: {'scores':{},'class4':{}}}
    missedFiles = 0
    
    # Load observations index and filter by location, time, and exclude list
    listCasts = pkg_obs.observations_load_and_filter(opt, "CTD", ana_start, ana_end)
    totalFiles = len(listCasts)
    if totalFiles == 0:
        print('No casts to process')
        return

    with open(logfile,'w') as f:
        f.write('Error messages from CTD cast analysis run at ' + str(datetime.datetime.utcnow()) + ':\n\n')

        print("Begin calculating metrics for {} casts ...".format(totalFiles))
        if opt['parallel']:
            with mp.Pool(processes=min(opt['nproc'], len(listCasts))) as pool:
                # flagging change - GO 20230612
                iterator = [pool.apply_async(process_a_cast_GO, args=(opt, cast, depth_levels, sdomains)) for cast in listCasts]
                results = [r.get() for r in iterator]
        else:
           # flagging change - GO 20230612
            results = [process_a_cast_GO(opt, cast, depth_levels, sdomains) for cast in listCasts]

        # Put the results in the output structure
        for result in results:
            scores, class4s, name, e = result
            if e is not None:
                if e != 'single pt':
                    missedFiles += 1
                    f.write(str(missedFiles) + '. Couldn\'t process cast ' + name + ' for ' + casename + ' because:\n        ' + str(e) + '\n\n')
            else:
                metrics[casename]['scores'][name],metrics[casename]['class4'][name] = scores, class4s

    print('Finished processing {} casts.'.format(totalFiles))
    if missedFiles > 0:
        print('{} casts could not be analyzed, see error log: {}'.format(missedFiles,logfile))

    return metrics

def process_a_cast(opt, cast, depth_levels):
    try:
        castPath = cast['filename']
        name = cast['code']
        obsdata = pkg_obs.load_observation(castPath, 'CTD')
        moddata, pfile = pkg_data.load_extraction(opt, cast)
        if moddata is None:
            # This cast must have been skipped during extraction
            return None,None, name, None
        scores, class4s, status = getSkill(obsdata, moddata, depth_levels)
        if status == 1:
            return scores, class4s, name, None
        else:
            return None, None, name, 'single pt'
    except Exception as e:
        return None, None, name, e

# added code to get subdomain - GO 20230612
def process_a_cast_GO(opt, cast, depth_levels, sdomains):

    try:
        
        castPath = cast['filename']
        name = cast['code']
        obsdata = pkg_obs.load_observation(castPath, 'CTD')
        pTemp_o,pSal_o,oDep = obsdata['pTemp'], obsdata['salinity'], obsdata['z']
        time = obsdata['time']
        pos = obsdata['lat'], obsdata['lon']
    
        # return none if obs data empty - GO 20230612
        if not oDep.any():
            #print("no observation founds in obs data file ", castPath)
            return None,None, name, None
        
        moddata, pfile = pkg_data.load_extraction(opt, cast)
        if moddata is None:
            # This cast must have been skipped during extraction
            return None,None, name, None

        # returns additional metrics (mean, std) in scores and includes subdomain - GO 20230612
        scores, class4s, status = getSkill_GO(obsdata, moddata, depth_levels, sdomains)
            
        if status == 1:
            return scores, class4s, name, None
        else:
            return None, None, name, 'single pt'
    except Exception as e:
        print("exception: ", e)
        return None, None, name, e
