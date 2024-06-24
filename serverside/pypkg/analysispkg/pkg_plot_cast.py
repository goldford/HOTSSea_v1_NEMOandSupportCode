import cmocean.cm as cm
import datetime
import gsw
import logging
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import pickle
import warnings
import shapely.geometry

from analysispkg import pkg_analyze_cast
from analysispkg import pkg_utils
from analysispkg import pkg_plot

plt.switch_backend('agg')
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.colors as colors
from matplotlib.path import Path as Path
import matplotlib.gridspec as gridspec

# GO added 20230509
#import matplotlib.cm as cm
import matplotlib as mpl
#from scipy.ndimage.filters import gaussian_filter

# added some plot variables for fine tuning plot params -GO
# to override / activate my changes - GO 2023-06-05
go_changes = 'on'


def plot_casts(opt):
    paths = {'analysis_base': os.path.join(opt['dir_process'], 'CTD')}
    minimum_samples_per_region = opt['plot']['CTD']['minimum_samples_per_region']
    formats = opt['plot']['file_formats']

    # Get bounding boxes and map zooms
    print('Loading domain file and computing bounding boxes...')
    depth_levels, coords = pkg_analyze_cast.read_CTD_domain(opt['analysis']['CTD']['domain_file'])
    bb = utils().get_bounding_boxes(coords)

    bathy = pkg_plot.get_bathy(opt)

    for period, (ana_start, ana_end) in opt['analysis']['periods'].items():
        print("Working on period {} ({} to {})".format(period, ana_start, ana_end))

        paths['figure_output_base_folder'] = os.path.join(opt['dir_plots'], 'CTD', period)
        os.makedirs(paths['figure_output_base_folder'], exist_ok=True)

        analysis_times = {'start': ana_start, 'end': ana_end}

        pfile = os.path.join(paths['analysis_base'], 'CTDcast_metrics_{}.pickle'.format(period))

        # load previously saved validation metrics (from process_cast_analyze)
        print('Loading validation metrics from ' + pfile + '...')
        try:
            metrics = {}
            with open(pfile, 'rb') as f:
                add_metrics = pickle.load(f, encoding="bytes")
            if add_metrics is None:
                print('No metrics recorded in the file, skipping plotting for the period')
                continue
            else:
                for k in add_metrics.keys():
                    metrics[k] = add_metrics[k]
        except Exception as e:
            raise ValueError(f'Could not load validation data from {pfile}, because: {e}. '
                             f'Either fix the path or re-run analyze.py')

        # skipping b/c broken 20230616
        if go_changes != 'on':# make key plot for analysis
          print('Creating key plot...')
          plotters().makeKeyPlot(paths, bathy, bb, metrics, formats)
          print(' ')

        if go_changes == 'on':
         # plot temperature and salinity deviations as probability density
         print('Plotting profile PDF and means...')
         with warnings.catch_warnings():
             warnings.filterwarnings(action='ignore', message='Mean of empty slice')
             plotters().plotProfilePDFsAndMeans(paths, metrics, bb, minimum_samples_per_region, formats)
             # plot PDFs with outliers excluded
             std_thres = opt['analysis']['CTD'].get('std_thres', None)
             if std_thres is not None:
                 plotters().plotProfilePDFsAndMeans(paths, metrics, bb, minimum_samples_per_region, formats, std_thres=std_thres)
         print(' ')
        
        else:
          # unmodified by GO
         # plot temperature and salinity deviations as probability density
         print('Plotting profile PDFs...')
         with warnings.catch_warnings():
             warnings.filterwarnings(action='ignore', message='Mean of empty slice')
             plotters().plotProfilePDFs(paths, metrics, bb, minimum_samples_per_region, formats)
             # plot PDFs with outliers excluded
             std_thres = opt['analysis']['CTD'].get('std_thres', None)
             if std_thres is not None:
                 plotters().plotProfilePDFs(paths, metrics, bb, minimum_samples_per_region, formats, std_thres=std_thres)
         print(' ')
          
        
        # summary maps which include RMSE by depth class
        if go_changes != 'on':
            # Not working for ORAS5 - omitting GO 20230418
            # plot summary maps of regional bias and RMS error
            print('Plotting summary maps...')
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                warnings.filterwarnings(action='ignore',
                                        message='Data has no positive values, and therefore cannot be log-scaled.')
                plotters().plotSummaryMaps(opt, paths, metrics, bathy, bb, minimum_samples_per_region, analysis_times)
            print(' ')

            # OMITTING - GO 20230418
            # plot all individual casts
            print('Plotting individual casts. This might take a while...')
            plotters().plotAllCasts(opt, paths, metrics, bb, bathy, depth_levels)

    print('Finished.')


def plot_casts_comparisons(opt, mods):
    paths = {}
    minimum_samples_per_region = opt['CTD']['minimum_samples_per_region']
    formats = opt['file_formats']

    # Get bounding boxes
    print('Loading domain file and computing bounding boxes...')
    depth_levels, coords = pkg_analyze_cast.read_CTD_domain(opt['CTD']['domain_file'])
    bb = utils().get_bounding_boxes(coords)

    bathy = pkg_plot.get_bathy(opt)

    print('Generating overview map ...')
    ## fac controls the buffer zone around the edges of the polygon, i.e. fac = 0. means polygon edges are on edges of plot
    # _, map_overview = make_map_overview(coords, fac=0.01)

    for period in opt['analysis']['periods'].keys():
        print("Working on period {}".format(period))

        # add metadata from runs that we are comparing to
        meta = {}
        picklefiles = []
        runKeys = []
        summaryKeys = []
        summaryColors = []
        summaryMarkers = []

        paths['figure_output_base_folder'] = os.path.join(opt['dir_plots'], 'CTD', period)
        os.makedirs(paths['figure_output_base_folder'], exist_ok=True)

        #//logfile = os.path.join(paths['analysis_base'], 'CTDcast_analysis_errorLog_{}.txt'.format(period))

        for mi,mod in enumerate(mods):
            runid = mod['runid']

            meta[runid] = mod
            picklefile = 'CTDcast_metrics_%s.pickle'%period
            picklefiles += [os.path.join(mod['dir_process'], 'CTD', picklefile)]

            runKeys.append(mod['casename'])
            summaryKeys.append(mod['casename'])
            summaryColors.append(opt['mods'][mi]['color'])
            summaryMarkers.append(opt['mods'][mi]['marker'])

        try:
            metrics = {}
            for picklefile in picklefiles:
                print('Loading data from %s...'%picklefile)
                with open(picklefile,'rb') as f:
                    add_metrics = pickle.load(f,encoding="bytes")
                    if add_metrics is not None:
                        for k in add_metrics.keys():
                            metrics[k] = add_metrics[k]
                    else:
                        print('No data from %s, skipping...'%runid)
        except Exception as e:
            raise ValueError('Could not load validation data from %s, because: %s. Either fix the path or re-run analyze.py'%(picklefile,e))

        # plot comparisons of temperature and salinity deviations as probability density for runs specified by summaryKeys
        print('Plotting comparisons of profile PDFs...')
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore', message='Mean of empty slice')
            plotters().plotProfilePDF_comparisons(paths, opt, metrics, bb, minimum_samples_per_region, summaryKeys,
                                                  summaryColors, period, formats)
        print(' ')

        # #plot TS distance statistics
        # print('Plotting TS distance statistics...')
        # plotters().plot_TS_distance(paths,metrics,coords,runKeys,summaryColors,summaryMarkers)
        # print(' ')

        # make key plot for analysis
        print('Creating key plot...')
        plotters().makeKeyPlot_comparisons(paths, bathy, bb, metrics, formats, opt)
        print(' ')

    print('Finished.')


@pkg_utils.time_it
def make_map_zooms(opt,coords,fac = 0.01):
    zooms = sorted(coords.keys())
    if opt['parallel']:
        # run in parallel
        with mp.Pool(processes=min(opt['nproc'], len(zooms)+1)) as pool:
            # Make the per-zoom maps
            iterator = [pool.apply_async(make_map_zoom, args=(coords, zoom, fac)) for zoom in zooms]
            # Make the overview map
            iterator += [pool.apply_async(make_map_overview, args=(coords, fac))]
            results = [r.get() for r in iterator]
    else:
        results = [make_map_zoom(coords, zoom, fac) for zoom in zooms]
        results += [make_map_overview(coords, fac)]

    maps = {zoom: map for zoom,map in results}

    return maps


def make_map_overview(coords, fac):
    lats = np.array([])
    lons = np.array([])
    for zoom in sorted(coords.keys()):
        lats = np.hstack([lats,coords[zoom][:,0]])
        lons = np.hstack([lons,coords[zoom][:,1]])
    map = Basemap(projection='merc',
                  llcrnrlat=(1. - fac) * lats.min(),
                  llcrnrlon=(1. - np.sign(lons.min()) * fac) * lons.min(),
                  urcrnrlat=(1. + fac) * lats.max(),
                  urcrnrlon=(1. + np.sign(lons.max()) * fac) * lons.max(),
                  resolution='h')
    return "overview", map


def make_map_zoom(coords, zoom, fac):
    latmin = (1. - fac) * coords[zoom][:,0].min()
    latmax = (1. + fac) * coords[zoom][:,0].max()
    lonmin = (1. - np.sign(coords[zoom][:,1].min()) * fac) * coords[zoom][:,1].min()
    lonmax = (1. + np.sign(coords[zoom][:,1].min()) * fac) * coords[zoom][:,1].max()
    d = np.sqrt((latmax - latmin)**2. + (lonmax - lonmin)**2.)
    if d < 2.:
        res = 'f'
    elif d < 10.:
        res = 'h'
    else:
        res = 'i'
    map = Basemap(projection = 'merc',llcrnrlat = latmin,llcrnrlon = lonmin,
                  urcrnrlat = latmax,urcrnrlon = lonmax,resolution = res)
    return zoom, map


class utils():
    def get_bounding_boxes(self,coords):
        bb = {}
        for k in coords.keys():
            # make sure the polygon is closed
            if np.all(coords[k][-1] == coords[k][0]):
                bb[k] = coords[k]
            else:
                bb[k] = np.vstack([coords[k], coords[k][0, :]])
        return bb

    def getMonRange(self,yr,analysis_times):
        if yr == analysis_times['start'].year:
            if yr != analysis_times['end'].year:
                return range(analysis_times['start'].month,13)
            else:
                return range(analysis_times['start'].month,analysis_times['end'].month + 1)
        elif yr == analysis_times['end'].year:
            return range(1,analysis_times['end'].month + 1)
        else:
            return range(1,13)

    def getYearMonKeys(self,yr,mon,metrics):
        out = []
        for f in metrics['scores'].keys():
            #TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
            if metrics['scores'][f] is None:
                continue
            mDT = metrics['scores'][f]['ObsTime'][0]
            if np.logical_and(int(mDT.year) == yr,int(mDT.month) == mon):
                out.append(f)

        return out

    def getAvg_BiasRMSE(self,metrics,bb,yr,mon,dKey,var):
        path = Path(bb)
        keys = utils().getYearMonKeys(yr,mon,metrics)
        if len(keys) == 0:
            return np.nan,np.nan,0

        scores = metrics['scores']

        write = np.zeros([4,])
        for f in keys:
            # TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
            if scores[f] is None:
                continue
            pos = scores[f]['ObsLocation']
            bias = scores[f]['bias'][dKey][var]
            rmse = scores[f]['rmse'][dKey][var]
            write = np.vstack([write,np.hstack([pos,bias,rmse])])

        write = write[1:,:]

        c = np.where(path.contains_points(write[:,:2]))[0]
        write = write[c,:]

        N = len(c)
        bias = np.nanmean(write[:,2])
        rmse = np.nanmean(write[:,3])

        return bias,rmse,N

    def get_common_profiles(self, metrics, opt):
        profiles = []
        for model_name, m in metrics.items():
            profiles.append(np.array([[profile_name, c4['dep'].max()]
                                      for profile_name, c4 in m['class4'].items()
                                      if c4 is not None]).astype(str))

        # get common profile names with intersect-union approach
        plist = [p[:, 0] for p in profiles]
        common_profiles = pkg_utils.get_common_stations(opt, metrics.keys(), plist)
        common_profiles = np.array(list(common_profiles))

        # get minimum common bottom depth for each profile
        common_depths = np.zeros([len(common_profiles), ]).astype(float)
        for i, profile in enumerate(common_profiles):
            depths = []
            for p in profiles:
                ip = p[:, 0] == profile
                if np.any(ip):
                    depths.append(float(p[ip, 1]))
            common_depths[i] = np.min(depths)

        return common_profiles, common_depths

    # added by GO 20230607, started w/ copy of get_profile_pdfs
    # notes - not handling outliers, sample size filtering rn (ignore std_thres, min_sample)
    def get_profile_means(self, scores, class4, bb, min_sample, reg, common_profiles, common_depths, std_thres=None):
        """

        Parameters
        ----------
        scores : dict
        class4 : dict
        bb : dict
            CTD region boundaries
        min_sample : int
            Minimum number of samples per layer. Exclude layers with fewer samples.
        reg : str
            CTD region name (from the CTD_analysis_domain...yaml)
        std_thres : float, optional
            Threshold for filtering out the outlier profiles. Outliers are detected as model-obs diff values
            outside of `mean +- thres*stdev` range. No outlier detection by default.

        Returns
        -------
        
        {'mean_obs': [], 'mean_mod': [], 'std_obs': [], 'std_mod': []}

        """
        bounds = Path(bb[reg])
        # profiles within the region
        # c4b = [c4 for k, c4 in class4.items() if bounds.contains_point(scores[k]['ObsLocation'])]
        # TODO (2020-10-01): Check why Nones appear, fix, go back to the above
        #kc4b = [(k, c4) for k, c4 in class4.items()
        #        if scores[k] is not None and bounds.contains_point(scores[k]['ObsLocation'])]
        if common_profiles is None:
            kc4b = [(k, c4) for k, c4 in class4.items()
                    if scores[k] is not None and bounds.contains_point(scores[k]['ObsLocation'])]
        else:
            kc4b = [(k, c4) for k, c4 in class4.items()
                    if str(k) in common_profiles and scores[k] is not None
                    and bounds.contains_point(scores[k]['ObsLocation'])]
                    
        print("n profiles: ", len(kc4b))
        # early exit if no valid profiles
        if len(kc4b) == 0:
            #return {'dep': [], 'bins': [], 'pdf': [], 'n': [], 'raw': [], 'outliers': []}, \
            #       {'dep': [], 'bins': [], 'pdf': [], 'n': [], 'raw': [], 'outliers': []}
            return {'mean_obs': [], 'mean_mod': [], 'std_obs': [], 'std_mod': []}, \
                   {'mean_obs': [], 'mean_mod': [], 'std_obs': [], 'std_mod': []}
        kb = [_[0] for _ in kc4b]  # profile names
        c4b = [_[1] for _ in kc4b]  # class4 values

        if common_profiles is None:
            z = [c4['dep'] for c4 in c4b]  # depth vectors for each profile
        else:
            z = [c4['dep'][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
        lens = [len(_) for _ in z]  # profile lengths
        maxlen = max(lens)
        z = z[np.argmax(lens)]
        # tails of nans for each profile to make them same length
        tails = [np.full(maxlen - len1, np.nan) for len1 in lens]
        
        def var_prep(var): 
        
            """ Get diffs and determine outliers """
            if common_profiles is None:
                
                # added by GO
                mvs = [c4['model'][var] for c4 in c4b]
                ovs = [c4['obs'][var] for c4 in c4b]
            
                #v = [c4['model'][var] - c4['obs'][var] for c4 in c4b]  # collect model-obs diffs
            else:
                mvs = [c4['model'][var][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
                ovs = [c4['obs'][var][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
                
                #v = [mv - ov for mv,ov in zip(mvs,ovs)]

            #v = [np.append(_, tail) for _, tail in zip(v, tails)]  # extend with nans to unify lengths
            #v = np.column_stack(v)  # array (ndep, nprofiles)
            #v_abs = np.abs(v)
            
            # added by GO
            mvs = [np.append(_, tail) for _, tail in zip(mvs, tails)]  # extend with nans to unify lengths
            mvs = np.column_stack(mvs)  # array (ndep, nprofiles)
            ovs = [np.append(_, tail) for _, tail in zip(ovs, tails)]  # extend with nans to unify lengths
            ovs = np.column_stack(ovs)  # array (ndep, nprofiles)
            
            
            # not sure how to translate over yet - GO
            # (filters if too few samples found at each depth)      
            #binlim = np.ceil(np.nanmax(v_abs))  # bounds for the histogram
            #nv = np.sum(v_abs < binlim, axis=1)  # nans not counted here

            # exclude undersampled layers
            #z_well = nv >= min_sample
            #v = v[z_well, :]
            #zv = z[z_well]
                        
            # not handling outliers here - GO
            # TODO Return a message
            # determine outliers
            #if (std_thres is None) | (std_thres == "None"):
            #    iout = np.full(v.shape[1], False)
            #    out_names = []
            #else:
            
            #    v_std = np.nanstd(v, axis=1, keepdims=True)
            #    v_mean = np.nanmean(v, axis=1, keepdims=True)
            #    
            #    iout = np.abs(v - v_mean) > std_thres * v_std  # NOTE nans are not counted as outliers here
            #    iout = np.any(iout, axis=0)  # profiles with outliers at ANY layer are marked as "bad"
            #    out_names = [kb1 for kb1, i1 in zip(kb, iout) if i1]
            
            mvs_mean = np.nanmean(mvs, axis=1, keepdims=True)
            mvs_std = np.nanstd(mvs, axis=1, keepdims=True)
            
            ovs_mean = np.nanmean(ovs, axis=1, keepdims=True)
            ovs_std = np.nanstd(ovs, axis=1, keepdims=True)
                
            #return v, zv, ~iout, out_names
            return ovs_mean, mvs_mean, ovs_std, mvs_std
            
            

        def var_pdf(v, zv):
            """ PDF (normalized histogram) along dim=1 of a 2D array """
            if v.size == 0:  # shortcut if number of profiles is below min_sample
                return [], [], [], [], []
            binlim = np.ceil(np.nanmax(np.abs(v)))  # bounds for the histogram
            
            # GO -20230607
            num_bins = 100
            
            if go_changes == 'on':
              if num_bins > v.shape[1]:
                num_bins = v.shape[1]
              bins = np.linspace(-1.0 * binlim, binlim, num_bins)  # histogram bins
            else:
              bins = np.linspace(-1.0 * binlim, binlim, v.shape[1])  # histogram bins
            
            pdf = [np.histogram(v1, bins=bins)[0] for v1 in v]  # histogram for each depth
            pdf = np.vstack(pdf)  # (ndep, nbins)
            nv = pdf.sum(axis=1)  # hist sum for each depth
            pdf = pdf / nv[:, None]  # normalized histogram
            # inz = nv != 0  # do not divide by zero
            # pdf[inz, :] = pdf[inz, :] / nv[inz][:, None]  # normalized histogram

            # exclude undersampled layers
            z_good = nv >= min_sample
            if not z_good.any():
                zv = []
                # v = np.zeros((0, 0)) + np.nan
            else:
                nv = nv[z_good]
                zv = zv[z_good]
                pdf = pdf[z_good, :]
                v = v[z_good, :]
                # v[v == 0.] = np.nan  # why this?
            return zv, bins, pdf, nv, v

        
        #t, zt, igood_t, out_t = var_prep('T') 
        #s, zs, igood_s, out_s = var_prep('S')
        
        t_obs_mean, t_mod_mean, t_obs_std, t_mod_std = var_prep('T') 
        s_obs_mean, s_mod_mean, s_obs_std, s_mod_std = var_prep('S') 
        
        #igood_ts = igood_t & igood_s
        #zt, bins_t, pdf_t, nt, t = var_pdf(t[:, igood_ts], zt)
        #zs, bins_s, pdf_s, ns, s = var_pdf(s[:, igood_ts], zs)

        return {'mean_obs': t_obs_mean, 'mean_mod': t_mod_mean, 'std_obs': t_obs_std, 'std_mod': t_mod_std}, \
               {'mean_obs': s_obs_mean, 'mean_mod': s_mod_mean, 'std_obs': s_obs_std, 'std_mod': s_mod_std}
               
    def get_profile_pdfs(self, scores, class4, bb, min_sample, reg, common_profiles, common_depths, std_thres=None):
        """

        Parameters
        ----------
        scores : dict
        class4 : dict
        bb : dict
            CTD region boundaries
        min_sample : int
            Minimum number of samples per layer. Exclude layers with fewer samples.
        reg : str
            CTD region name (from the CTD_analysis_domain...yaml)
        std_thres : float, optional
            Threshold for filtering out the outlier profiles. Outliers are detected as model-obs diff values
            outside of `mean +- thres*stdev` range. No outlier detection by default.

        Returns
        -------

        """
        bounds = Path(bb[reg])
        # profiles within the region
        # c4b = [c4 for k, c4 in class4.items() if bounds.contains_point(scores[k]['ObsLocation'])]
        # TODO (2020-10-01): Check why Nones appear, fix, go back to the above
        #kc4b = [(k, c4) for k, c4 in class4.items()
        #        if scores[k] is not None and bounds.contains_point(scores[k]['ObsLocation'])]
        if common_profiles is None:
            kc4b = [(k, c4) for k, c4 in class4.items()
                    if scores[k] is not None and bounds.contains_point(scores[k]['ObsLocation'])]
        else:
            kc4b = [(k, c4) for k, c4 in class4.items()
                    if str(k) in common_profiles and scores[k] is not None
                    and bounds.contains_point(scores[k]['ObsLocation'])]
        # early exit if no valid profiles
        if len(kc4b) == 0:
            return {'dep': [], 'bins': [], 'pdf': [], 'n': [], 'raw': [], 'outliers': []}, \
                   {'dep': [], 'bins': [], 'pdf': [], 'n': [], 'raw': [], 'outliers': []}
        kb = [_[0] for _ in kc4b]  # profile names
        c4b = [_[1] for _ in kc4b]  # class4 values

        if common_profiles is None:
            z = [c4['dep'] for c4 in c4b]  # depth vectors for each profile
        else:
            z = [c4['dep'][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
        lens = [len(_) for _ in z]  # profile lengths
        maxlen = max(lens)
        z = z[np.argmax(lens)]
        # tails of nans for each profile to make them same length
        tails = [np.full(maxlen - len1, np.nan) for len1 in lens]

        def var_prep(var):
            """ Get diffs and determine outliers """
            if common_profiles is None:
                v = [c4['model'][var] - c4['obs'][var] for c4 in c4b]  # collect model-obs diffs
            else:
                mvs = [c4['model'][var][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
                ovs = [c4['obs'][var][c4['dep'] <= common_depths[common_profiles == profile_name][0]] for profile_name,c4 in zip(kb,c4b)]
                v = [mv - ov for mv,ov in zip(mvs,ovs)]

            v = [np.append(_, tail) for _, tail in zip(v, tails)]  # extend with nans to unify lengths
            v = np.column_stack(v)  # array (ndep, nprofiles)
            v_abs = np.abs(v)
            binlim = np.ceil(np.nanmax(v_abs))  # bounds for the histogram
            nv = np.sum(v_abs < binlim, axis=1)  # nans not counted here
            # exclude undersampled layers
            z_well = nv >= min_sample
            v = v[z_well, :]
            zv = z[z_well]
            # TODO Return a message
            # determine outliers
            if (std_thres is None) | (std_thres == "None"):
                iout = np.full(v.shape[1], False)
                out_names = []
            else:
                v_std = np.nanstd(v, axis=1, keepdims=True)
                v_mean = np.nanmean(v, axis=1, keepdims=True)
                iout = np.abs(v - v_mean) > std_thres * v_std  # NOTE nans are not counted as outliers here
                iout = np.any(iout, axis=0)  # profiles with outliers at ANY layer are marked as "bad"
                out_names = [kb1 for kb1, i1 in zip(kb, iout) if i1]
            return v, zv, ~iout, out_names

        def var_pdf(v, zv):
            """ PDF (normalized histogram) along dim=1 of a 2D array """
            if v.size == 0:  # shortcut if number of profiles is below min_sample
                return [], [], [], [], []
            binlim = np.ceil(np.nanmax(np.abs(v)))  # bounds for the histogram
            
            # GO -20230607
            num_bins = 100
            
            if go_changes == 'on':
              if num_bins > v.shape[1]:
                num_bins = v.shape[1]
              bins = np.linspace(-1.0 * binlim, binlim, num_bins)  # histogram bins
            else:
              bins = np.linspace(-1.0 * binlim, binlim, v.shape[1])  # histogram bins
            
            pdf = [np.histogram(v1, bins=bins)[0] for v1 in v]  # histogram for each depth
            pdf = np.vstack(pdf)  # (ndep, nbins)
            nv = pdf.sum(axis=1)  # hist sum for each depth
            pdf = pdf / nv[:, None]  # normalized histogram
            # inz = nv != 0  # do not divide by zero
            # pdf[inz, :] = pdf[inz, :] / nv[inz][:, None]  # normalized histogram

            # exclude undersampled layers
            z_good = nv >= min_sample
            if not z_good.any():
                zv = []
                # v = np.zeros((0, 0)) + np.nan
            else:
                nv = nv[z_good]
                zv = zv[z_good]
                pdf = pdf[z_good, :]
                v = v[z_good, :]
                # v[v == 0.] = np.nan  # why this?
            return zv, bins, pdf, nv, v

        t, zt, igood_t, out_t = var_prep('T')
        s, zs, igood_s, out_s = var_prep('S')
        igood_ts = igood_t & igood_s
        zt, bins_t, pdf_t, nt, t = var_pdf(t[:, igood_ts], zt)
        zs, bins_s, pdf_s, ns, s = var_pdf(s[:, igood_ts], zs)

        return {'dep': zt, 'bins': bins_t, 'pdf': pdf_t, 'n': nt, 'raw': t, 'outliers': out_t}, \
               {'dep': zs, 'bins': bins_s, 'pdf': pdf_s, 'n': ns, 'raw': s, 'outliers': out_s}

    def makeProfilePDFs(self, scores, class4, bb, minSample, reg):
        # find max diff for T and S and depth for all profiles
        maxDiffT = 0.
        maxDiffS = 0.
        maxDep = np.array([])
        count = 0
        for k in class4.keys():
            # TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
            if scores[k] is None:
                continue
            # pick only belonging to this region
            if not Path(bb[reg]).contains_point(scores[k]['ObsLocation']):
                continue
            diffT = np.abs(class4[k]['model']['T'] - class4[k]['obs']['T'])
            if np.nanmax(diffT) > maxDiffT:
                maxDiffT = np.nanmax(diffT)

            diffS = np.abs(class4[k]['model']['S'] - class4[k]['obs']['S'])
            if np.nanmax(diffS) > maxDiffS:
                maxDiffS = np.nanmax(diffS)

            if len(class4[k]['dep']) > len(maxDep):
                maxDep = class4[k]['dep']
            count += 1

        diffGridT = np.linspace(-1. * np.ceil(maxDiffT), np.ceil(maxDiffT), count)
        xT, yT = np.meshgrid(diffGridT, maxDep)
        diffGridS = np.linspace(-1. * np.ceil(maxDiffS), np.ceil(maxDiffS), count)
        xS, yS = np.meshgrid(diffGridS, maxDep)

        # calculate PDFs (normalized histograms) for T and S diffs
        pdf_T = np.zeros_like(xT)
        pdf_S = np.zeros_like(xS)
        count = 0
        for k in class4.keys():  # go over profiles
            # TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
            if scores[k] is None:
                continue
            if not Path(bb[reg]).contains_point(scores[k]['ObsLocation']):
                continue

            diffT = class4[k]['model']['T'] - class4[k]['obs']['T']
            T1 = np.zeros_like(maxDep)
            T1[:len(diffT)] = diffT
            diffS = class4[k]['model']['S'] - class4[k]['obs']['S']
            S1 = np.zeros_like(maxDep)
            S1[:len(diffS)] = diffS
            if count == 0:
                T = T1
                S = S1
            else:
                T = np.column_stack([T,T1])
                S = np.column_stack([S,S1])

            for i in range(len(diffT)):  # for each depth level
                # TODO Bug: This counts nans and places them in the last bin
                indT = np.digitize(diffT[i], diffGridT) - 1
                # if i == len(maxDep)-1 and indT == len(diffGridT)-1:
                #     print(f"=====count: {count}    diffT[i]: {diffT[i]}")
                indS = np.digitize(diffS[i], diffGridS) - 1
                pdf_T[i, indT] += 1
                pdf_S[i, indS] += 1
            count += 1

        N = np.sum(pdf_T, axis=1)
        for i in range(pdf_T.shape[0]):
            pdf_T[i,:] = pdf_T[i,:]/np.sum(pdf_T[i,:])
            pdf_S[i,:] = pdf_S[i,:]/np.sum(pdf_S[i,:])

        pdf_T = pdf_T[N >= minSample,:]
        pdf_S = pdf_S[N >= minSample,:]
        maxDep = maxDep[N >= minSample]
        xT,yT = np.meshgrid(diffGridT,maxDep)
        xS,yS = np.meshgrid(diffGridS,maxDep)
        if not (N >= minSample).any():
            T = np.zeros((0,0)) + np.nan
            S = np.zeros((0,0)) + np.nan
        else:
            T = T[N >= minSample,:]
            S = S[N >= minSample,:]
            N = N[N >= minSample]

            T[T == 0.] = np.nan
            S[S == 0.] = np.nan

        return {'x': xT, 'y': yT, 'pdf': pdf_T, 'raw': T}, {'x': xS, 'y': yS, 'pdf': pdf_S, 'raw': S}, N, maxDep

    def plotRhoConts(self,ax):
        sGrid = np.linspace(0.,36.,1000)
        tGrid = np.linspace(-4.,30.,1000)
        s,t = np.meshgrid(sGrid,tGrid)
        pGrid = gsw.rho(s,t,np.zeros_like(s))
        c = ax.contour(sGrid,tGrid,pGrid - 1000.,np.arange(0.,30.,2.),colors='k')
        plt.clabel(c,fmt='%.1f',inline=1)

    def getAnnualTS(self,ts):
        arr = {}
        for src in ['obs','model']:
            arr[src] = np.zeros([2,])
            for mon in ts.keys():
                arr[src] = np.vstack([arr[src],ts[mon][src]])
            arr[src] = arr[src][1:,:]

        return arr

    def getAxLims(self,ts):
        lims = {}
        for d in ['x','y']:
            lims[d] = {}
            lims[d]['min'] = 10.**10
            lims[d]['max'] = -10.**10
        arr = np.zeros([2,])
        for run in sorted(ts.keys()):
            ts2 = utils().getAnnualTS(ts[run])
            for src in ['obs','model']:
                ts2[src][ts2[src] > 100.] = np.nan
            for i,d in enumerate(['y','x']):
                lims[d]['min'] = np.min(np.array([lims[d]['min'],np.nanmin(ts2['obs'][:,i]),np.nanmin(ts2['model'][:,i])]))
                lims[d]['max'] = np.max(np.array([lims[d]['max'],np.nanmax(ts2['obs'][:,i]),np.nanmax(ts2['model'][:,i])]))

        lims['x']['min'] = 0.95 * lims['x']['min']
        if lims['y']['min'] > 0.:
            lims['y']['min'] = 0.95 * lims['y']['min']
        else:
            lims['y']['min'] = 1.05 * lims['y']['min']
        lims['x']['max'] = 1.05 * lims['x']['max']
        lims['y']['max'] = 1.05 * lims['y']['max']

        return lims

    def make_TS_dict(self,data,coords,runKeys):
        TS = {}
        for reg in sorted(coords.keys()):
            first_cast_flag = True
            for cast in data[runKeys[0]]['scores'].keys():
                for run in runKeys:
                    if cast not in data[run]['class4'].keys():
                        continue
                    elif Path(coords[reg]).contains_point(data[run]['scores'][cast]['ObsLocation']):
                        if first_cast_flag:
                            TS[reg] = {}
                            first_cast_flag = False
                        if run not in TS[reg].keys():
                            TS[reg][run] = {}
                        date = '%s/%s'%(str(data[run]['scores'][cast]['ObsTime'][0].year),\
                                        str(data[run]['scores'][cast]['ObsTime'][0].month).zfill(2))
                        if date not in TS[reg][run].keys():
                            TS[reg][run][date] = {}
                            start = {'obs':True,'model':True}
                        for src in ['obs','model']:
                            write = np.column_stack([data[run]['class4'][cast][src]['T'],\
                                                     data[run]['class4'][cast][src]['S']])
                            if start[src]:
                                TS[reg][run][date][src] = write
                                start[src] = False
                            else:
                                TS[reg][run][date][src] = np.vstack([TS[reg][run][date][src],write])

        return TS

    def get_monrange_TS(self,TS):
        dates = []
        for reg in sorted(TS.keys()):
            for run in TS[reg].keys():
                dates.extend([m for m in TS[reg][run].keys()])
        dates = sorted(dates)
        minmon = dates[0]; mindate = datetime.datetime(int(minmon[:minmon.find('/')]),int(minmon[minmon.find('/')+1:]),15)
        maxmon = dates[-1]; maxdate = datetime.datetime(int(maxmon[:maxmon.find('/')]),int(maxmon[maxmon.find('/')+1:]),15)

        printmons = ['%s/%s'%(str(mindate.year),str(mindate.month).zfill(2))]
        date = mindate
        while date < maxdate:
            yr = date.year
            mon = date.month
            if mon == 12:
                mon = 1
                yr = yr + 1
            else:
                mon = mon + 1
            date = datetime.datetime(yr,mon,15)
            printmons.append('%s/%s'%(str(date.year),str(date.month).zfill(2)))

        return printmons

    def get_n_subplots(self,n,base):
        if n == 1:
            return 1, 1, (base, base)

        j = int(np.ceil(np.sqrt(n)))
        i = int(np.floor(j%n)) - 1
        if i*j < n:
            i += 1
        i = np.max([i,1])

        aspect_ratio = float(j)/float(i)
        j_size = aspect_ratio * base
        i_size = base / aspect_ratio
        size = (i_size,j_size)

        return j,i,size


def plot_single_cast(configs, class4s, bathy, locs, f, cols, depth_levels, outdir, formats):
    fig = plt.figure(figsize=(12, 6), dpi=200, constrained_layout=True)
    subfigs = fig.subfigures(1, 2, width_ratios=[1, 2])  # , wspace=None, hspace=None

    # Plot the map
    axm = subfigs[0].subplots(2, 1)
    pkg_plot.plot_bathy_with_zoom(axm, bathy, locs[f])
    # show transect direction: arrow from 1st to last point
    x, y = locs[f]['lon'], locs[f]['lat']
    axm[1].scatter(x, y, facecolor='r', edgecolor='k', lw=2, zorder=1001)
    # axm[1].text(0.98, 0.98, f"{f}\n{locs[f]['time']}", ha='right', va='top', transform=axm[1].transAxes)
    axm[1].text(0.98, 0.98, locs[f]['time'], ha='right', va='top', transform=axm[1].transAxes)  # cast time
    # axm[1].text(0.98, 0.98, f, ha='right', va='top', transform=axm[1].transAxes)  # cast name, includes time

    (axT, axS) = subfigs[1].subplots(1, 2, sharey=True)
    axT.set_ylabel('Pressure (dBar)')

    alph = 0.5
    for ax, var, xlab in zip([axT, axS], ['T', 'S'], [r'Potential Temperature ($^\circ$C)', 'Salinity (psu)']):
        handles = []
        labels = []
        for a, config, class4, col in zip(range(len(cols)), sorted(configs), class4s, cols):
            if class4 is None:
                print('Couldn\'t find ' + var + ' in ' + f + ' for ' + config)
                continue

            ax.plot(class4['model'][var], class4['dep'], color=col, lw=3, alpha=alph)
            if a == 0:
                ax.plot(class4['obs'][var], class4['dep'], color='k', lw=3, ls='--')
                handles.append(Line2D([], [], color='k', lw=3, ls='--'))
                labels.append('Observation')
            if var == 'T':
                handles.append(Line2D([], [], color=col, lw=3, alpha=alph))
                labels.append(config)

        ax.grid()
        ax.set_xlabel(xlab)

        if var == 'T':
            ax.legend(handles, labels, loc=4)

        for dep in depth_levels:
            if dep[1] < class4['dep'].max():
                ax.axhline(dep[1], color='k', ls='--')

    axT.invert_yaxis()  # also inverts for axS as they are linked

    svPath = os.path.join(outdir, f + '_class4')
    pkg_plot.figure_saver(svPath, formats, fig)
    print('Saved ' + f + ' to ' + svPath)
    plt.close()


def save_outliers(out_t, out_s, save_path):
    """ Print and save a list of profile names in the lists """
    def print_list(names, var):
        if len(names) > 0:
            print(f'----- Outliers for {var}:')
            for outn in names:
                print(outn)
            print('----- end of outliers')

    out_ts = list(set(out_t + out_s))  # merge lists
    print_list(out_t, 'temperature')
    print_list(out_s, 'salinity')
    print_list(out_ts, 'temperature or salinity')
    # save a yaml-ready list in a file
    if len(out_ts) > 0:
        with open(save_path, 'w') as fid:
            fid.write('[')
            for out in out_ts[:-1]:
                fid.write(f"'*/{out}*'," + '\n')
            fid.write(f"'*/{out_ts[-1]}*']" + '\n')
    # TODO Are profile names guaranteed to be unique?


class plotters():
    def makeKeyPlot(self, paths, bathy, bb, metrics, formats):
        keys = [k for k in metrics.keys()]
        dummy_key = keys[0]
        posx, posy = [], []
        for K in sorted(metrics[dummy_key]['scores'].keys()):
            # TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
            if metrics[dummy_key]['scores'][K] is not None:
                pos = metrics[dummy_key]['scores'][K]['ObsLocation']
                posx += [pos[1]]
                posy += [pos[0]]
        self.plotKeyPlot(paths, bathy, bb, posx, posy, formats)

    def makeKeyPlot_comparisons(self, paths, bathy, bb, metrics, formats, opt):
        # get location of all casts
        locations = {}
        for _, cmetric in metrics.items():
            for stn, stnobj in cmetric['scores'].items():
                if stnobj is not None:
                    locations[stn] = stnobj['ObsLocation']
        # extract just the locations from the common profiles
        common_profiles, _ = utils().get_common_profiles(metrics, opt)
        posx, posy = [], []
        for stn in common_profiles:
            pos = locations[stn]
            posx += [pos[1]]
            posy += [pos[0]]
        self.plotKeyPlot(paths, bathy, bb, posx, posy, formats)

    def plotKeyPlot(self, paths, bathy, bb, posx, posy, formats):
        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        _ = pkg_plot.plot_bathy(ax, bathy, label_axes=True)
        
        if go_changes != 'on':

            ax.scatter(posx, posy, s=20, facecolor='r', edgecolor='k', lw=1, zorder=1000)
            for k in sorted(bb.keys()):
              x, y = bb[k][:, 1], bb[k][:, 0]
              ax.plot(x, y, color='k', zorder=1001)
              # area centroid
              centr = shapely.geometry.Polygon(bb[k][:, ::-1]).centroid
              xa, ya = centr.x, centr.y
              # ax.scatter(xa, ya, color='r', zorder=1002)  # large red dot for each area
              ht = ax.text(xa, ya, k, zorder=1003, weight='bold', ha='center', va='bottom', color='k', fontsize=8)  # opaque text
              #ht.set_bbox(dict(facecolor='w', alpha=0.5))  # half-transparent white background around text
        else:
        
            #GO 20230509 commented out
            #ht.set_bbox(dict(facecolor='w', alpha=0.5))  # half-transparent white background around text
            
            #ax.scatter(posx, posy, s=5, facecolor='r', edgecolor='r', lw=0, zorder=1000, alpha=0.3) #GO 20230509 
        
            #GO 20230509 - does not quite work
#           heatmap, xedges, yedges = np.histogram2d(posx, posy, bins=1000)
#           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]
#           heatmap = gaussian_filter(heatmap, sigma=16)
#           ax.imshow(heatmap,extent=extent, cmap=cm.amp,  zorder=1000)

            #GO 20230509 - hexbin attempt
            # https://stackoverflow.com/questions/48613920/use-of-extend-in-a-pcolormesh-plot-with-discrete-colorbar
            # discrete colors - create x boundaries to have x-1 intervals
            boundaries = np.array([1, 25, 50, 100, 200, 400, 800, 1200])
            # create list of colors from colormap
            cmap_amp = plt.cm.get_cmap('cmo.amp',len(boundaries))
            colors = list(cmap_amp(np.arange(len(boundaries))))
            #replace first color with white
            #colors[0] = "white"
            cmap = mpl.colors.ListedColormap(colors[:-1], "")
            # set over-color to last color of list 
            cmap.set_over(colors[-1])
            norm = mpl.colors.BoundaryNorm(boundaries, ncolors=len(boundaries)-1, clip=False)
        
            hb = ax.hexbin(posx, posy, gridsize=100, cmap=cmap, norm=norm, mincnt=1, extent=[min(posx), max(posx), min(posy), max(posy)])
            cb = fig.colorbar(hb, ax=ax)
      
            for k in sorted(bb.keys()):
              x, y = bb[k][:, 1], bb[k][:, 0]
              ax.plot(x, y, color='k', zorder=1001)
              # area centroid
              centr = shapely.geometry.Polygon(bb[k][:, ::-1]).centroid
              xa, ya = centr.x, centr.y
              # ax.scatter(xa, ya, color='r', zorder=1002)  # large red dot for each area
              ht = ax.text(xa, ya, k, zorder=1003, weight='bold', ha='center', va='bottom', color='k', fontsize=8)  # opaque text
            
              #GO 20230509 commented out
              #ht.set_bbox(dict(facecolor='w', alpha=0.5))  # half-transparent white background around text

        save_path = os.path.join(paths['figure_output_base_folder'], 'keyPlot')
        msg = pkg_plot.figure_saver(save_path, formats, fig, dpi=200, bbox_inches='tight')
        print(msg)

    def plotAllCasts(self, opt, paths, metrics, bb, bathy, depth_levels):
        months = ['', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
                  'November', 'December']
        formats = opt['plot']['file_formats']
        cols = [np.array([1, 0, 0])]

        if opt['parallel']:
            pool = mp.Pool(processes=opt['nproc'])
        iterator = []

        dummyKey = sorted(metrics.keys())[0]
        for reg in sorted(bb.keys()):
            bb_Path = Path(bb[reg])
            # m = maps[reg]

            for mon in range(1,13):
                outdir = os.path.join(paths['figure_output_base_folder'], 'Profiles', reg, months[mon])
                os.makedirs(outdir, exist_ok=True)

                keys = []
                locs = {}
                for f in sorted(metrics[dummyKey]['class4'].keys()):
                    # TODO (2020-10-01): Check why Nones appear, fix, remove the if-check
                    if metrics[dummyKey]['scores'][f] is None:
                        continue
                    if metrics[dummyKey]['scores'][f]['ObsTime'][0].month == mon:
                        if bb_Path.contains_point(metrics[dummyKey]['scores'][f]['ObsLocation']):
                            locs[f] = {}
                            locs[f]['lat'] = metrics[dummyKey]['scores'][f]['ObsLocation'][0]
                            locs[f]['lon'] = metrics[dummyKey]['scores'][f]['ObsLocation'][1]
                            locs[f]['time'] = metrics[dummyKey]['scores'][f]['ObsTime'][0]
                            keys.append(f)

                if len(keys) == 0:
                    print('No data for ' + months[mon] + ' in ' + reg)
                    continue
                else:
                    for f in keys:
                        class4s = []
                        for config in metrics.keys():
                            try:
                                class4 = metrics[config]['class4'][f]
                            except:
                                class4 = None
                            class4s += [class4]

                        configs = list(metrics.keys())
                        if opt['parallel']:
                            iterator += [pool.apply_async(plot_single_cast,
                                                          args=(configs, class4s, bathy, locs, f, cols, depth_levels,
                                                                outdir, formats))]
                        else:
                            plot_single_cast(configs, class4s, bathy, locs, f, cols, depth_levels, outdir, formats)

        if opt['parallel']:
            _ = [r.get() for r in iterator]
            pool.close()
            pool.join()

    def plotSummaryMaps(self, opt, paths, metrics, m, bb, minSample, analysis_times):
        logging.basicConfig(level=logging.INFO)

        formats = opt['plot']['file_formats']

        # generate random color palette for bar graph showing sample numbers
        # can be replaced by specifying cols manually as a list of color definitions
        i = 0
        cols = []
        while i < len(bb.keys()):
            cols.append(np.random.rand(3,))
            i += 1

        if opt['parallel']:
            pool = mp.Pool(processes=opt['nproc'])
        iterator = []

        # loop through available model runs
        for mKey in sorted(metrics.keys()):
            print('Processing data for ' + mKey + '...')

            # loop through depth levels
            keys = [k for k in metrics[mKey]['scores'].keys()]
            if len(keys) == 0:
                print("No metrics for {}, moving on".format(mKey))
                continue
            dummyKey = keys[0]
            
            # none catching bit - GO 20230410    
            go1 = 0
            if metrics[mKey]['scores'][dummyKey] is None:
                print(metrics[mKey]['scores'])
                while metrics[mKey]['scores'][dummyKey] is None:
                    go1 += 1
                    if go1 < len(keys):
                        dummyKey = keys[go1]
                    else:
                        print("cannot find score for model / instrument")
                        break
        
            for dKey in sorted(metrics[mKey]['scores'][dummyKey]['bias'].keys()):
                if dKey == 'Full':
                    depStr = 'fullDepth'
                elif dKey[0] == '>':
                    depStr = 'below' + dKey[1:] + 'm'
                else:
                    depStr = dKey.replace('->', 'to') + 'm'

                outdirt = os.path.join(paths['figure_output_base_folder'], 'Summaries', 'Temperature', depStr)
                outdirs = os.path.join(paths['figure_output_base_folder'], 'Summaries', 'Salinity', depStr)
                os.makedirs(outdirt, exist_ok=True)
                os.makedirs(outdirs, exist_ok=True)

                for yr in range(analysis_times['start'].year, analysis_times['end'].year + 1):
                    monRange = utils().getMonRange(yr, analysis_times)
                    for mon in monRange:
                        for var, longvar, units, outdir in \
                                zip(['T', 'S'], ['Temperature', 'Salinity'], [r'$^\circ$C', 'psu'], [outdirt, outdirs]):
                            args = metrics[mKey], m, bb, yr, mon, mKey, dKey, var, longvar, \
                                   units, minSample, cols, outdir, depStr, formats
                            if opt['parallel']:
                                iterator += [pool.apply_async(plotters().plotSummaryMapHelper, args=args)]
                            else:
                                plotters().plotSummaryMapHelper(*args)
        if opt['parallel']:
            _ = [r.get() for r in iterator]
            pool.close()
            pool.join()

    def plotSummaryMapHelper(self, metric, bathy, bb, yr, mon, mKey, dKey, var, longvar, units, minSample, cols, outdir,
                             depStr, formats):
        try:
            months = ['', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

            baseSize = 440
            rmse_amp = 1.  # scaling parameters for size of dots
            bias_amp = 3.  # bias colorbar scale

            fig = plt.figure(figsize=(6,7))

            m_ax = plt.subplot(1,1,1)
            N_ax = inset_axes(m_ax,width='30%',height='20%',loc=1)

            plt.sca(m_ax)
            _ = pkg_plot.plot_bathy(m_ax, bathy, label_axes=True)

            locs = np.zeros([2,])
            stats = np.zeros([2,])
            labels = sorted(bb.keys())
            for i,reg in enumerate(labels):
                locs = np.vstack([locs,bb[reg].mean(axis=0)])
                
                # GO 20230411 - left here, should allow for entire hindcast analysis, entire year not just each year, month
                bias,rmse,N = utils().getAvg_BiasRMSE(metric,bb[reg],yr,mon,dKey,var)
                
                
                if N < minSample:
                    bias = np.nan
                    rmse = np.nan
                stats = np.vstack([stats,np.hstack([bias,rmse])])
                N_ax.bar(i,N,0.8,align='center',facecolor=cols[i],edgecolor='k',lw=2)
            locs = locs[1:,:]
            stats = stats[1:,:]

            x, y = locs[:, 1], locs[:, 0]

            pc = m_ax.scatter(x, y, c=stats[:, 0], edgecolor='k', lw=2, s=(rmse_amp * stats[:, 1] * baseSize),
                              zorder=1001, cmap=cm.balance, vmax=bias_amp, vmin=-1. * bias_amp)
            m_ax.text(0.02, 0.02, longvar + '\n' + months[mon] + '/%1i' % yr, ha='left', va='bottom',
                      transform=m_ax.transAxes)
            plt.colorbar(pc, ax=m_ax, shrink=0.5, pad=0.05,
                         label='E[' + var + r'$_{model}$ - ' + var + r'$_{obs}$] (' + units + ')')

            hand = m_ax.scatter(x[-1], y[-1], facecolor='0.5', edgecolor='k', s=(rmse_amp * 0.25 * baseSize), lw=2,
                                zorder=1000, label='RMSE = %.1f' % (rmse_amp * 0.25) + ' ' + units)
            m_ax.legend(loc=2)
            hand.set_visible(False)

            N_ax.set_xticks(np.arange(0., len(bb.keys()), 1))
            N_ax.set_xticklabels(labels, rotation=90.)
            N_ax.set_yscale('log')
            N_ax.set_ylabel('# of casts')
            N_ax.set_ylim(1., 200.)
            N_ax.set_yticks([1., 10., 100.])
            for l in N_ax.get_yticks():
                N_ax.axhline(l, color='k', ls='--', lw=0.5)

            save_path = os.path.join(outdir, f"{var}_Summary_{depStr}_{months[mon]}{str(yr)}_{mKey}")
            pkg_plot.figure_saver(save_path, formats, fig, dpi=300, bbox_inches='tight')
            plt.close()

            logging.info(f"Saved {save_path}")

        except Exception as e:
            logging.error(e, exc_info=True)


    # added by GO 20230607
    def plotProfilePDFsAndMeans(self, paths, metrics, bb, min_sample, formats, std_thres=None):
            
      
        # manual setting of x-axis limits
        xLimT = 3.
        xLimS = 5.       #Go 20230411 commented line to use these below back in
        lw1 = 1        #GO 20230411 CRMSE
        lw2 = 1        #GO 20230411  +/- lines
        fs1 = 8          #GO 20230411  fontsize legend
        gamma1 = 0.65      #Go 20230411 for powernorm color norm
        fig_dim_1 = 9
        fig_dim_2 = 3
        min_x_salt_mean = 22 # for plotting, found through trial and error across all plots
        max_x_salt_mean = 34
        min_x_temp_mean = 6
        max_x_temp_mean = 13
        grd_hspace = 0.4 # percent of axis vert reserved for labels
        grd_wspace = 0.3
        

          
        outdir = os.path.join(paths['figure_output_base_folder'], 'ProfilePDFsAndMeans')
        os.makedirs(outdir, exist_ok=True)

        for run in metrics.keys():
            print('Processing data from ' + run + '...')
            for reg in bb.keys():
            
                print('Plotting ' + reg + '...')
                
                # returns pcolormeshes for plots
                plotT,plotS,N,maxDep = utils().makeProfilePDFs(metrics[run]['scores'],metrics[run]['class4'],bb,min_sample,reg)
                
                
                plotT, plotS = utils().get_profile_pdfs(metrics[run]['scores'], metrics[run]['class4'], bb, min_sample, reg, None, None, std_thres)
                
                if len(plotT['dep']) == 0 and len(plotS['dep']) == 0:
                    print(f'Region {reg} has no profiles. Not making a plot for this')
                    continue
                
                
                # GO added below
                plotT_mean, plotS_mean = utils().get_profile_means(metrics[run]['scores'], metrics[run]['class4'], bb, min_sample, reg, None, None, std_thres)
                # for now take dep from fnc above
                plotT_mean['dep'] = plotT['dep']
                plotS_mean['dep'] = plotS['dep']
                
                # temporarily turn this off - GO 20230607
                #save_outliers(plotT['outliers'], plotS['outliers'], os.path.join(outdir, f"{reg}_outliers.txt"))
                
                fig = plt.figure(figsize=(fig_dim_1, fig_dim_2))

                #gs = gridspec.GridSpec(1, 6)
                gs = gridspec.GridSpec(8, 10,hspace=grd_hspace,wspace=grd_wspace)
                
                
                axt_mean = plt.subplot(gs[0:5,:2])
                axt = plt.subplot(gs[0:5,2:4],sharey=axt_mean)
                axs_mean = plt.subplot(gs[0:5,4:6], sharey=axt_mean)
                axs = plt.subplot(gs[0:5,6:8], sharey=axt_mean)
                
                ax = [axt_mean, axt, axs_mean, axs]   
                
                # plot the mean temps
                ax[0].set_ylabel('Pressure (dBar)', fontsize=fs1)
                dep_n = len(plotT_mean['dep'])
                obs_t_mean = plotT_mean['mean_obs'][0:dep_n]
                ax[0].plot(obs_t_mean, plotT_mean['dep'], 'k-')
                obs_t_mean = plotT_mean['mean_mod'][0:dep_n]
                ax[0].plot(obs_t_mean, plotT_mean['dep'], 'k--')
                ax[0].set_xlabel("Mean\nTemp. ($^\circ$C)", fontsize=fs1)
                
                mean_obs = Line2D([], [], color='k', lw=lw2)
                mean_mod = Line2D([], [], color='k', ls='--', lw=lw2)
                ax[0].legend([mean_obs, mean_mod], ['Obs.', 'Mod.'], fontsize=fs1, bbox_to_anchor=[0.4, -0.4])
                #ax[0].legend([mean_obs, mean_mod], ['Obs.', 'Mod.'], fontsize=fs1, loc='upper left')
                
                
                # plot the mean salin
                dep_n = len(plotS_mean['dep'])
                obs_s_mean = plotS_mean['mean_obs'][0:dep_n]
                ax[2].plot(obs_s_mean, plotS_mean['dep'], 'k-')
                obs_s_mean = plotS_mean['mean_mod'][0:dep_n]
                ax[2].plot(obs_s_mean, plotS_mean['dep'], 'k--')
                ax[2].set_xlabel("Mean\nSalin. (PSU)", fontsize=fs1)
                
                
                # don't plot values deeper than 2000 m, arbitrarily based on depth of ARGO profiles
                zMax = np.min(np.array([2000., max(plotT['dep'][-1], plotS['dep'][-1])]))

                xlabs = [u'$\Delta$ Temp.\nMod. - Obs. ($^\circ$C)',
                         u'$\Delta$ Salin.\nMod. - Obs. (psu)']
                
                #set up the pdfs:
                for plot, a, xLim, xLab in zip([plotT, plotS], [ax[1],ax[3]], [xLimT, xLimS], xlabs):
                    if len(plot['dep']) == 0:
                        continue
                    # pc = a.pcolor(plot['x'], plot['y'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                    #               cmap=cm.amp, norm=colors.LogNorm(vmin=0.01, vmax=1.), zorder=1000, shading='flat')
                    bin_centers = plot['bins'][:-1] + np.diff(plot['bins']) / 2
                    
                    
                    #pc = a.pcolor(bin_centers, plot['dep'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                    #              cmap=cm.amp, norm=colors.PowerNorm(gamma=gamma1), shading='nearest',rasterized=True)
                    
                    pc = a.pcolor(bin_centers, plot['dep'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                                  cmap=cm.amp, norm=colors.LogNorm(vmin=0.01, vmax=1.), shading='nearest',
                                  rasterized=True)
                    
                    
                    vmean = np.nanmean(plot['raw'], axis=1) # can the util code above to be modified to do real vmean? -GO
                    vstd = np.nanstd(plot['raw'], axis=1)
                    a.plot(vmean, plot['dep'], color='k', lw=lw1, label='Bias')
                    a.plot(vmean - 1. * vstd, plot['dep'], color='k', ls='--', lw=lw2, label='CRMSE')
                    a.plot(vmean + 1. * vstd, plot['dep'], color='k', ls='--', lw=lw2)
                    a.set_ylim(0., zMax)
                    a.grid(which='both')
                    a.set_xlim(-1.*xLim, xLim)
                    a.set_xlabel(xLab, fontsize=fs1)
                    
                    
                bias = Line2D([], [], color='k', lw=lw2)
                rms = Line2D([], [], color='k', ls='--', lw=lw2)
                ax[1].legend([bias, rms], ['Bias', 'CRMSE'], fontsize=fs1, bbox_to_anchor=[0.2, -0.4])
                #ax[1].legend([bias, rms], ['Bias', 'CRMSE'], fontsize=fs1, loc='upper left')
                
                # # of obs plot
                ax.append(plt.subplot(gs[0:5,8], sharey=axt_mean))
                # not sure what this is for -GO
                #ax.append(fig.add_axes([0.8, 0.3, 0.02, 0.4]))
                ax[4].plot(plotT['n'], plotT['dep'], color='k', label='T')
                
                if len(plotT['n']) != len(plotS['n']) or np.any(plotT['n'] != plotS['n']):
                    ax[4].plot(plotS['n'], plotS['dep'], '--', color='k', label='S')
                    #ax[4].legend(fontsize=fs1)
                ax[4].set_xlabel('# of \nobservations', fontsize=fs1)
                ax[4].grid(which='both')
                ax[4].set_ylim(0., zMax)
                ax[4].set_xlim(0, 1.05 * max(plotT['n'].max(), plotS['n'].max()))

                ax[0].invert_yaxis()

                # colorbar
                ax_clr = plt.subplot(gs[-1,3:9])
                cbar = plt.colorbar(pc, cax=ax_clr, shrink=0.3, orientation='horizontal', pad=2)
                
#                cbar = plt.colorbar(pc, cax=ax[-1], label='Relative Density of Observations')                
#                #cbar = plt.colorbar(pc, label='Relative density of obs.', ax=ax[0:2], shrink=0.6)
                cbar.set_ticks([0.01, 0.05, 0.1, 0.25, 0.5, 1.])
                cbar.set_ticklabels([r'$\leq$1%', '5%', '10%', '25%', '50%', '100%'])
                cbar.ax.tick_params(labelsize=fs1)
                cbar.ax.set_xlabel('Relative density of observations', fontdict={'fontsize':fs1})

                # title w/ region and model run
                #ax[4].text(1.18, 0.98, reg + '\n' + run, ha='left', va='top', transform=ax[2].transAxes)
                # text above is in data units for position so use title
                
                # since there's space above plots, use it:
                if reg == "DI":
                  subtitle = "Discovery Is. (DI)" 
                elif reg == "SGN":
                  subtitle = "Strait of Georgia North (SGN)"
                elif reg == "SGS":
                  subtitle = "Strait of Georgia South (SGS)"
                elif reg == "GI":
                  subtitle = "Southern Gulf Islands (GI)"
                elif reg == "HS":
                  subtitle = "Haro Strait (HS)"
                elif reg == "PS":
                  subtitle = "Puget Sound (PS)"
                elif reg == "JFS":
                  subtitle = "Juan de Fuca Strait (JFS)"
                else:
                  subtitle = reg
                
                ax[2].set_title(subtitle,fontdict={'fontsize':fs1},loc='left')

                # format the pdf's
#                for a in [ax[0], ax[1]]:
#                    a.axvline(0., color='k', ls='dotted', lw=lw2, zorder=1003)
#                    lims = a.get_xlim()
#                    start = np.floor(lims[1]/0.5)*0.5
#                    ticks = np.linspace(-1. * start, start, 7)
#                    print(ticks)
#                    a.set_xticks(ticks)
#                    plt.sca(a)
#                    plt.tick_params(axis='x', labelsize=fs1)
#                    plt.tick_params(axis='y', labelsize=fs1)
                
                # format axis for temp pdf
                ax[1].axvline(0., color='k', ls='dotted', lw=lw2, zorder=1003)
                lims = ax[1].get_xlim()
                start = np.floor(lims[1]/0.5)*0.5
                ticks = np.linspace(-1. * start, start, 7)
                ax[1].set_xticks(ticks)
                plt.sca(ax[1])
                plt.tick_params(axis='x', labelsize=fs1)
                plt.tick_params(axis='y', labelsize=fs1)
                
                # format axis for salt pdf
                ax[3].axvline(0., color='k', ls='dotted', lw=lw2, zorder=1003)
                lims = ax[3].get_xlim()
                #start = np.floor(lims[1]/0.5)*0.5
                #ticks = np.linspace(-1. * start, start, 7)
                ticks = [-4,-2,0,2,4] # overriding salt
                
                ax[3].set_xticks(ticks)
                plt.sca(ax[3])
                plt.tick_params(axis='x', labelsize=fs1)
                plt.tick_params(axis='y', labelsize=fs1)
                
                # format the axes for means
                ax_n = 0
                for a in [ax[0], ax[2]]:
                    
                    
                    if ax_n == 0:
                      a.set_xlim(min_x_temp_mean,max_x_temp_mean)
                      ticks = np.arange(min_x_temp_mean, max_x_temp_mean+1, 2)
                    else:
                      a.set_xlim(min_x_salt_mean,max_x_salt_mean)
                      ticks = np.arange(min_x_salt_mean, max_x_salt_mean+3, 3)
                    a.set_xticks(ticks)
                    a.grid(which='both')
                    plt.sca(a)
                    plt.tick_params(axis='x', labelsize=fs1)
                    plt.tick_params(axis='y', labelsize=fs1)
                    ax_n+=1
                    
                for a in [ax[1], ax[2], ax[3], ax[4]]:
                    plt.sca(a)
                    plt.tick_params(axis='y', labelleft=False)
                    plt.tick_params(axis='x', labelsize=fs1)

                # Tag the plots if outliers were excluded.
                # With no detected outliers, the plots from the "no-outlier" run are overwritten
                if len(plotT['outliers']) > 0 or len(plotS['outliers']) > 0:
                    tag = '_without_outliers'
                else:
                    tag = ''
                save_path = os.path.join(outdir, reg + f'_ProfilePDFAndMean{tag}')
                pkg_plot.figure_saver(save_path, formats, fig, dpi=450, bbox_inches='tight')
                plt.close()
      


    def plotProfilePDFs(self, paths, metrics, bb, min_sample, formats, std_thres=None):
        
        if go_changes == 'on':
          # manual setting of x-axis limits
          xLimT = 3.
          xLimS = 5.       #Go 20230411 commented line to use these below back in
          lw1 = 0.5        #GO 20230411 CRMSE
          lw2 = 0.5        #GO 20230411  +/- lines
          fs1 = 6          #GO 20230411  fontsize legend
          gamma1 = 0.65      #Go 20230411 for powernorm color norm
          fig_dim_1 = 8
          fig_dim_2 = 5
        else:
          # manual setting of x-axis limitsi've go
          lw1 = 3.
          lw2 = 3
          xLimT = 3.
          xLimS = 3.
          fs1 = 9

        
          
        outdir = os.path.join(paths['figure_output_base_folder'], 'ProfilePDFs')
        os.makedirs(outdir, exist_ok=True)

        for run in metrics.keys():
            print('Processing data from ' + run + '...')
            for reg in bb.keys():
            
                print('Plotting ' + reg + '...')
                
                plotT,plotS,N,maxDep = utils().makeProfilePDFs(metrics[run]['scores'],metrics[run]['class4'],bb,min_sample,reg)
                plotT, plotS = utils().get_profile_pdfs(metrics[run]['scores'], metrics[run]['class4'], bb, min_sample, reg, None, None, std_thres)
                if len(plotT['dep']) == 0 and len(plotS['dep']) == 0:
                    print(f'Region {reg} has no profiles. Not making a plot for this')
                    continue

                save_outliers(plotT['outliers'], plotS['outliers'], os.path.join(outdir, f"{reg}_outliers.txt"))
                
                if go_changes == 'on':
                  fig = plt.figure(figsize=(fig_dim_1, fig_dim_2))
                else:
                  fig = plt.figure(figsize=(8, 5))

                gs = gridspec.GridSpec(1, 6)
                axt = plt.subplot(gs[:2])
                axs = plt.subplot(gs[2:4], sharey=axt)
                ax = [axt, axs]
                # don't plot values deeper than 2000 m, arbitrarily based on depth of ARGO profiles
                zMax = np.min(np.array([2000., max(plotT['dep'][-1], plotS['dep'][-1])]))

                xlabs = [u'$\Delta$ Temperature\nModel - Observations ($^\circ$C)',
                         u'$\Delta$ Salinity\nModel - Observations (psu)']
                for plot, a, xLim, xLab in zip([plotT, plotS], ax, [xLimT, xLimS], xlabs):
                    if len(plot['dep']) == 0:
                        continue
                    # pc = a.pcolor(plot['x'], plot['y'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                    #               cmap=cm.amp, norm=colors.LogNorm(vmin=0.01, vmax=1.), zorder=1000, shading='flat')
                    bin_centers = plot['bins'][:-1] + np.diff(plot['bins']) / 2
                    
                    if go_changes == 'on':
                      pc = a.pcolor(bin_centers, plot['dep'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                                  cmap=cm.amp, norm=colors.PowerNorm(gamma=gamma1), shading='nearest',
                                  rasterized=True)
                    else:
                      pc = a.pcolor(bin_centers, plot['dep'], np.ma.masked_where(plot['pdf'] == 0., plot['pdf']),
                                  cmap=cm.amp, norm=colors.LogNorm(vmin=0.01, vmax=1.), shading='nearest',
                                  rasterized=True)
                    
                    
                    vmean = np.nanmean(plot['raw'], axis=1)
                    vstd = np.nanstd(plot['raw'], axis=1)
                    a.plot(vmean, plot['dep'], color='k', lw=lw1, label='Bias')
                    a.plot(vmean - 1. * vstd, plot['dep'], color='k', ls='--', lw=lw2, label='CRMSE')
                    a.plot(vmean + 1. * vstd, plot['dep'], color='k', ls='--', lw=lw2)
                    a.set_ylim(0., zMax)
                    a.grid(which='both')
                    a.set_xlim(-1.*xLim, xLim)
                    a.set_xlabel(xLab, fontsize=fs1)
                    

                ax[0].set_ylabel('Pressure (dBar)', fontsize=fs1)

                ax.append(plt.subplot(gs[4], sharey=axt))
                ax.append(fig.add_axes([0.8, 0.3, 0.02, 0.4]))

                ax[2].plot(plotT['n'], plotT['dep'], color='k', label='T')
                if len(plotT['n']) != len(plotS['n']) or np.any(plotT['n'] != plotS['n']):
                    ax[2].plot(plotS['n'], plotS['dep'], '--', color='k', label='S')
                    ax[2].legend(fontsize=fs1)
                    
                if go_changes == 'on':
                  ax[2].set_xlabel('# of observations', fontsize=fs1)
                else: 
                  ax[2].set_xlabel('# of observations')
                  
                ax[2].grid(which='both')
                ax[2].set_ylim(0., zMax)
                ax[2].set_xlim(0, 1.05 * max(plotT['n'].max(), plotS['n'].max()))

                ax[0].invert_yaxis()

                cbar = plt.colorbar(pc, cax=ax[-1], label='Relative Density of Observations')
                cbar.set_ticks([0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.])
                cbar.set_ticklabels([r'$\leq$1 %', '5 %', '10 %', '25 %', '50 %', '75 %', '100 %'])
                cbar.ax.tick_params(labelsize=fs1)
                
                # bias = Line2D([], [], color='k', lw=3)
                # rms = Line2D([], [], color='k', ls='--', lw=3)
                # ax[2].legend([bias, rms], ['Bias', 'CRMSE'], fontsize=8, bbox_to_anchor=[1.18, 0.15])
                ax[0].legend(fontsize=fs1)

                ax[2].text(1.18, 0.98, reg + '\n' + run, ha='left', va='top', transform=ax[2].transAxes)

                for a in [ax[1], ax[2]]:
                    plt.sca(a)
                    plt.tick_params(axis='y', labelleft=False,labelsize=fs1)
                    plt.tick_params(axis='x', labelsize=fs1)

                for a in [ax[0], ax[1]]:
                    a.axvline(0., color='k', ls='dotted', lw=lw2, zorder=1003)
                    lims = a.get_xlim()
                    start = np.floor(lims[1]/0.5)*0.5
                    ticks = np.linspace(-1. * start, start, 5)
                    a.set_xticks(ticks)
                    plt.sca(a)
                    plt.tick_params(axis='x', labelsize=fs1)
                    plt.tick_params(axis='y', labelsize=fs1)
                    

                # Tag the plots if outliers were excluded.
                # With no detected outliers, the plots from the "no-outlier" run are overwritten
                if len(plotT['outliers']) > 0 or len(plotS['outliers']) > 0:
                    tag = '_without_outliers'
                else:
                    tag = ''
                save_path = os.path.join(outdir, reg + f'_ProfilePDF{tag}')
                pkg_plot.figure_saver(save_path, formats, fig, dpi=450, bbox_inches='tight')
                plt.close()

    def plotProfilePDF_comparisons(self, paths, opt, metrics, bb, min_sample, runs, cols, period, formats):
        
        # params i added - GO 20230719
        regions_compare = ['DI', 'SGS', 'JFS'] # subset of regions (bb.keys()) for comparison -GO 20230719
        lw_mean = 1 # for tuning lw bias of pdf plot - GO 20230719
        lw_grid = 0.5
        fig_width = 6
        fig_height = 5.5
        fig_c_fs1 = 6 # font size reg text, axes labs
        fig_c_fs2 = 7 # font size titles
        letters = ['(a)', '(d)', '(b)', '(e)', '(c)', '(f)']
        
        outdir = os.path.join(paths['figure_output_base_folder'], 'ProfilePDF_comparisons')
        os.makedirs(outdir, exist_ok=True)
        
        
        # GO changes 20230719:
        # - change to plot multiple regions on one grid of plots
        if go_changes == 'on':
          
          show_range = False # to show fill between ucl and lcl
          custom_run_lbl = True
          run_lbl_mapping = {'SS1500-RUN141':'v0.14','SS1500-RUN181':'v0.18','SS1500-RUN100':'v0.1','SS1500-RUN160':'v0.16', 'SS500-201905':'SalishSeaCast'}
          
          # check each region for data - GO 20230719
          # (filters out regions w/ no data to help prep dynamically sized plot)
          # to do: redundant- make not so
          regions_compare_2 = []
          for reg in bb.keys():  
            print(reg)
            if reg in regions_compare:
              print('Checking ' + reg + '...')
            
              # collect data for plotting
              data = []
              num_casts = []
              common_profiles,common_depths = utils().get_common_profiles(metrics, opt)
              for k, m in metrics.items():  # for each model
                t, s = utils().get_profile_pdfs(m['scores'], m['class4'], bb, min_sample, reg, common_profiles, common_depths)  # PDFs for T and S
                if (len(t['n']) == 0 or np.isnan(t['raw']).all()) and \
                   (len(s['n']) == 0 or np.isnan(s['raw']).all()):
                    data.append(None)
                else:
                    data.append((t, s))
                    num_casts.append(max((t['n'].max(), s['n'].max())))

              if all(v is None for v in data):
                print(f'No profiles in region {reg}. Omitting from plotting')
                continue
              else:
                print(f'Found data in region {reg}')
                regions_compare_2.append(reg)
                
          num_regs = len(regions_compare_2) # number of regions to plot - GO change 20230719
          
           
          fig = plt.figure(figsize=(fig_width, fig_height))
          #gs = gridspec.GridSpec(60, 40 * num_regs) # I think this exceeds max num cols -GO
          gs = gridspec.GridSpec(68, 9 * num_regs)
          axs = []
          r = 0 # index of reg added to plot
          for reg in bb.keys():
            
            print(reg)
          
            if reg in regions_compare_2:
              print('Plotting ' + reg + '...')
              
              # collect data for plotting
              data = []
              num_casts = []
              common_profiles,common_depths = utils().get_common_profiles(metrics, opt)
              for k, m in metrics.items():  # for each model
                t, s = utils().get_profile_pdfs(m['scores'], m['class4'], bb, min_sample, reg, common_profiles, common_depths)  # PDFs for T and S
                if (len(t['n']) == 0 or np.isnan(t['raw']).all()) and \
                 (len(s['n']) == 0 or np.isnan(s['raw']).all()):
                  data.append(None)
                else:
                  data.append((t, s))
                  num_casts.append(max((t['n'].max(), s['n'].max())))

              if all(v is None for v in data):
                print(f'No profiles in region {reg}. Omitting from plotting')
                continue
           
              print(t.keys())
           
              # gridspec columns setup
              print(r)
              col_strt = r * 9
              col_end = r * 9 + 7
              #col_strt_nobs = (r * 8) + 6 # for layout w/ num obs plot to right of each
              #col_end_nobs = (r * 8) + 7
              
              #if r == 0:
              ax0 = plt.subplot(gs[  :25, col_strt:col_end])
              axs.append(ax0) # T
              ax1 = plt.subplot(gs[32:57, col_strt:col_end])
              axs.append(ax1) # S
               
              #else:
                # alt way of doing it
                #ax0 = fig.add_subplot(gs[  :25, col_strt:col_end], sharey=ax0_col1) # T
                #ax1 = fig.add_subplot(gs[30:55, col_strt:col_end], sharey=ax1_col1) # S
               # ax0 = plt.subplot(gs[  :25, col_strt:col_end], sharey=axs[0])
               # axs.append(ax0) # T
               # ax1 = plt.subplot(gs[30:55, col_strt:col_end], sharey=axs[1])
               # axs.append(ax1) # S
                

              
              # num obs plot
              #ax2 = fig.add_subplot(gs[  :38, col_strt_nobs:col_end_nobs])
          
              handles = []
              labels = []
              for run, col, d in zip(runs, cols, data):  # for each model (cols = colors here)
                if d is None:
                  print(f'No profiles for run {run} in region {reg}. Omitting from plotting')
                  continue
                for vdict, a in zip(d, [ax0, ax1]):  # for T and S
                  vmean = np.nanmean(vdict['raw'], axis=1)
                  vstd = np.nanstd(vdict['raw'], axis=1)
                  if show_range:
                    a.fill_betweenx(vdict['dep'], vmean - 1. * vstd, vmean + 1. * vstd, color=col, alpha=0.3)
                  a.plot(vmean, vdict['dep'], color=col, lw=lw_mean)
                  a.tick_params(axis='both', which='both',labelsize=fig_c_fs1)

#                # number of profiles
#                ax2.plot(d[0]['n'], d[0]['dep'], color=col)
#                if len(d[0]['n']) != len(d[1]['n']) or np.any(d[0]['n'] != d[1]['n']):
#                  print('Different number of observations for T and S; showing number of obs. for T')
#                  print(f"Num profiles: region {reg}\t run {run}:\t T: {np.max(d[0]['n'])}\t S: {np.max(d[1]['n'])}")

                handles.append(Line2D([], [], color=col, lw=3))
                
                # abbreviate if SalishSea - GO 20230719
                substring = 'SalishSea'
                run_label = run.replace(substring, 'SS')
                
                if custom_run_lbl: 
                  run_label = run_lbl_mapping[run_label]
                
                labels.append(run_label)
                #labels.append(run)
            
              # grids and labels
              for a in [ax0, ax1]:
                xlim = np.abs(a.get_xlim()).max()
                a.set_xlim(-xlim, xlim)
                a.set_yscale('log')
                a.invert_yaxis()
                if r == 0:
                  a.set_ylabel('Pressure (dBar)',labelpad=2, fontdict={'fontsize':fig_c_fs1})
                else:
                  a.set_ylabel('')
                #a.grid(which='both')
                a.grid(which='major', lw=lw_grid)
                a.axvline(0., color='k', ls='dotted', lw=lw_grid, zorder=1003)
                a.axis('tight')
                
                

              ax0.set_xlabel(r'$\Delta$ Temp. ($^\circ$C)', labelpad=1, fontdict={'fontsize':fig_c_fs1}) 
              # to add line break put space and then '\n' or else TEX messes up
              ax1.set_xlabel(r'$\Delta$ Salin. (PSU)', labelpad=1, fontdict={'fontsize':fig_c_fs1})

#              ax2.set_xlabel('# of observations')
#              ax2.grid(which='both')
#              ax2.invert_yaxis()
#              ax2.set_xlim(0, max(num_casts) * 1.05)
              
              # common legend for all panels
              #handles.append(Line2D([], [], color='k', lw=3))  # bias
              #labels.append('Bias')
              if show_range:
                handles.append(Patch(facecolor='k', edgecolor='none', alpha=0.3))  # rms
                labels.append(r'Bias $\pm$ CRMSE')
              
              if r == 0:
                #ax2.legend(handles,labels,fontsize=8,bbox_to_anchor=[1.315, 0.35])
                ax2 = fig.add_subplot(gs[62:, 0:]) # third row - GO
                ax2.legend(handles,labels,fontsize=fig_c_fs1, ncol=3, loc='center')
              
              # ax1.legend(handles,labels,fontsize=8,bbox_to_anchor=[1.315, 0.35])
              #ax1.legend(handles, labels, fontsize=8, bbox_to_anchor=[1.45, 0.25])
              
              #ax2.text(-0.4, -0.12, f'Region: {reg}\nPeriod: {period}', ha='left', va='top', transform=ax2.transAxes,
              #       fontsize=14)
              
              reg_lab = 'Region'
              if reg == 'DI': 
                reg_lab = 'Discovery Is.'
              elif reg == 'SGN': 
                reg_lab = 'St. of Georgia - North'
              elif reg == 'SGS':
                reg_lab = 'St. of Georgia - South'
              elif reg == 'GI':
                reg_lab = 'Gulf Is.'
              elif reg == 'HS':
                reg_lab = 'Haro St.'
              elif reg == 'PS':
                reg_lab = 'Puget Sound'
              elif reg == 'JFS':
                reg_lab = 'Juan de Fuca St.'
                
              #ax0.set_title(f'{reg_lab} ({reg})',fontdict={'fontsize':fig_c_fs2},loc='center')
              ax0.set_title(f'{reg}',fontdict={'fontsize':fig_c_fs2},loc='center', y=0.97)
              # added by GO 20240220
              ax1.set_title(f'{reg}',fontdict={'fontsize':fig_c_fs2},loc='center', y=0.97)


              # GO added 20240220
              ax0.text(0.1,1.03,letters[r],transform=ax0.transAxes, ha='center', color='k', fontdict={'fontsize':fig_c_fs2})
              ax1.text(0.1,1.03,letters[r+3],transform=ax1.transAxes, ha='center', color='k', fontdict={'fontsize':fig_c_fs2})

              # to do: use code from notebook to add new function that returns metrics for this reg and period (grouped by depths)

#              import matplotlib.table as tbl
#              ax2.set_axis_off()
#              table = tbl.table(ax2, cellText = val3, rowLabels = val2, colLabels = val1, rowColours = ["palegreen"] * 10,
#           colColours =["palegreen"] * 10,
#           cellLoc ='center', 
#           loc ='upper left')
#         
#       ax2.add_table(table)
#
#              table = tbl.table(cellText=['','','','','','',''], colLabels=["bias"], rowLabels=labels,rowLoc="center")
#              ax2.axis("tight")
              ax2.axis("off")

              r += 1


          save_path = os.path.join(outdir, 'MultRegn_ProfilePDF_comparison_GO')
          pkg_plot.figure_saver(save_path, formats, fig, dpi=400, bbox_inches='tight')
          plt.close()
              
              
              
              

         
        
        
        
        # original code - no GO changes - 20230719
        else: # if go_changes != on
          for reg in bb.keys():  # go over regions
            print('Plotting ' + reg + '...')

            # collect data for plotting
            data = []
            num_casts = []
            common_profiles,common_depths = utils().get_common_profiles(metrics, opt)
            for k, m in metrics.items():  # for each model
                t, s = utils().get_profile_pdfs(m['scores'], m['class4'], bb, min_sample, reg, common_profiles, common_depths)  # PDFs for T and S
                if (len(t['n']) == 0 or np.isnan(t['raw']).all()) and \
                   (len(s['n']) == 0 or np.isnan(s['raw']).all()):
                    data.append(None)
                else:
                    data.append((t, s))
                    num_casts.append(max((t['n'].max(), s['n'].max())))

            if all(v is None for v in data):
                print(f'No profiles in region {reg}. Omitting from plotting')
                continue
            fig = plt.figure(figsize=(6.5, 9))

            gs = gridspec.GridSpec(60, 40)
            ax0 = plt.subplot(gs[  :25, :30])
            ax1 = plt.subplot(gs[30:55, :30])
            ax2 = plt.subplot(gs[  :38, 34:])
            handles = []
            labels = []
            for run, col, d in zip(runs, cols, data):  # for each model
                if d is None:
                    print(f'No profiles for run {run} in region {reg}. Omitting from plotting')
                    continue
                for vdict, a in zip(d, [ax0, ax1]):  # for T and S
                    vmean = np.nanmean(vdict['raw'], axis=1)
                    vstd = np.nanstd(vdict['raw'], axis=1)
                    a.fill_betweenx(vdict['dep'], vmean - 1. * vstd, vmean + 1. * vstd, color=col, alpha=0.3)
                    a.plot(vmean, vdict['dep'], color=col, lw=3)

                # number of profiles
                ax2.plot(d[0]['n'], d[0]['dep'], color=col)
                if len(d[0]['n']) != len(d[1]['n']) or np.any(d[0]['n'] != d[1]['n']):
                    print('Different number of observations for T and S; showing number of obs. for T')
                print(f"Num profiles: region {reg}\t run {run}:\t T: {np.max(d[0]['n'])}\t S: {np.max(d[1]['n'])}")

                handles.append(Line2D([], [], color=col, lw=3))
                labels.append(run)

            # grids and labels
            for a in [ax0, ax1]:
                xlim = np.abs(a.get_xlim()).max()
                a.set_xlim(-xlim, xlim)
                a.set_yscale('log')
                a.invert_yaxis()
                a.set_ylabel('Pressure (dBar)')
                a.grid(which='both')
                a.axvline(0., color='k', ls='dotted', lw=2, zorder=1003)

            ax0.set_xlabel(r'$\Delta$ Temperature, Model - Observations ($^\circ$C)')
            ax1.set_xlabel(r'$\Delta$ Salinity, Model - Observations (psu)')

            ax2.set_xlabel('# of observations')
            ax2.grid(which='both')
            ax2.invert_yaxis()
            ax2.set_xlim(0, max(num_casts) * 1.05)

            # common legend for all panels
            handles.append(Line2D([], [], color='k', lw=3))  # bias
            handles.append(Patch(facecolor='k', edgecolor='none', alpha=0.3))  # rms
            labels.append('Bias')
            labels.append(r'Bias $\pm$ CRMSE')
            # ax1.legend(handles,labels,fontsize=8,bbox_to_anchor=[1.315, 0.35])
            ax1.legend(handles, labels, fontsize=8, bbox_to_anchor=[1.45, 0.25])

            ax2.text(-0.4, -0.12, f'Region: {reg}\nPeriod: {period}', ha='left', va='top', transform=ax2.transAxes,
                     fontsize=14)

            save_path = os.path.join(outdir, reg + '_ProfilePDF_comparison')
            pkg_plot.figure_saver(save_path, formats, fig, dpi=300, bbox_inches='tight')
            plt.close()

    class TS_distance():
        def plot_monthly_TS_diagrams(self,basePath,TS,coords,runKeys,summaryColors):
            os.makedirs(basePath,exist_ok=True)

            mons = utils().get_monrange_TS(TS)

            for reg in sorted(coords.keys()):
                regPath = os.path.join(basePath,reg)
                os.makedirs(regPath,exist_ok=True)
                print('Plotting TS diagrams for %s. Saving to %s'%(reg,regPath))

                lims = utils().getAxLims(TS[reg])
                for run,col in zip(runKeys,summaryColors):
                    nj,ni,size = utils().get_n_subplots(len(mons),7.5)
                    fig,axes = plt.subplots(nj,ni,sharex='col',sharey='row',figsize=size)

                    plotLeg = True
                    for mon,ax in zip(mons,axes.flatten()[:len(mons)]):
                        if plotLeg:
                            obs = Line2D([],[],ls='none',marker='o',color='k')
                            mod = Line2D([],[],ls='none',marker='o',color=col)
                            ax.legend([obs,mod],['Observations','Model (%s)'%run],loc=1)
                            plotLeg = False
                        ax.text(0.02,0.98,mon,ha='left',va='top',transform=ax.transAxes)
                        if mon not in TS[reg][run].keys():
                            ax.text(0.98,0.02,'NO DATA',ha='right',va='bottom',transform=ax.transAxes)
                            continue
                        ts = TS[reg][run][mon]
                        utils().plotRhoConts(ax)
                        ax.scatter(ts['obs'][:,1],ts['obs'][:,0],color='k',alpha=0.05)
                        ax.scatter(ts['model'][:,1],ts['model'][:,0],color=col,alpha=0.05)

                    for ax in axes[-1,:]:
                        ax.set_xlim(lims['x']['min'],lims['x']['max'])
                        ax.set_xlabel('Practical Salinity (psu)')

                    for ax in axes[:,0]:
                        ax.set_ylim(lims['y']['min'],lims['y']['max'])
                        ax.set_ylabel(u'Potential\nTemperature\n($^\circ$C)')

                    svPath = os.path.join(regPath,'TSdiagram_monthly_%s.png'%run)
                    plt.savefig(svPath,dpi=300,bbox_inches='tight')
                    plt.close()

        def plot_annual_TS_diagrams(self,basePath,TS,coords,runKeys,summaryColors):
            os.makedirs(basePath,exist_ok=True)

            mons = utils().get_monrange_TS(TS)
            nj,ni,size = utils().get_n_subplots(len(runKeys),7.5)

            for reg in sorted(coords.keys()):
                print('Plotting diagram for %s.'%reg)
                lims = utils().getAxLims(TS[reg])

                fig,axes = plt.subplots(nj,ni,sharex='col',sharey='row',figsize=size)
                try:
                    b = axes[0,0]
                except:
                    axes = np.asarray(axes)
                    axes = axes.reshape(axes.size,1)
                axes[0,0].text(0.02,0.98,reg,ha='left',va='top',transform=axes[0,0].transAxes)
                for run,col,ax in zip(runKeys,summaryColors,axes.flatten()[:len(runKeys)]):
                    utils().plotRhoConts(ax)
                    ts = utils().getAnnualTS(TS[reg][run])
                    ax.scatter(ts['obs'][:,1],ts['obs'][:,0],color='k',alpha=0.025)
                    ax.scatter(ts['model'][:,1],ts['model'][:,0],color=col,alpha=0.025)

                    obs = Line2D([],[],ls='none',marker='o',color='k')
                    mod = Line2D([],[],ls='none',marker='o',color=col)
                    ax.legend([obs,mod],['Observations','Model (%s)'%run])

                for ax in axes[-1,:]:
                    ax.set_xlim(lims['x']['min'],lims['x']['max'])
                    ax.set_xlabel('Practical Salinity (psu)')

                for ax in axes[:,0]:
                    ax.set_ylim(lims['y']['min'],lims['y']['max'])
                    ax.set_ylabel(u'Potential Temperature ($^\circ$C)')

                plt.tight_layout()

                svPath = os.path.join(basePath,'TSdiagram_annual_%s.png'%(reg))
                plt.savefig(svPath,dpi=300,bbox_inches='tight')
                plt.close()

        def plot_regional_TS_distance(self,saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit'):
            mons = utils().get_monrange_TS(TS)
            regs = [k for k in sorted(coords.keys())]

            nj,ni,size = utils().get_n_subplots(len(regs),7.5)

            fig,axes = plt.subplots(nj,ni,sharex='col',sharey='row',figsize=size)
            try:
                b = axes[0,0]
            except:
                axes = np.asarray(axes)
                axes = axes.reshape(axes.size,1)
            maxD = 0.; minD = 0.
            for ax,reg in zip(axes.flatten()[:len(regs)],regs):
                ax.set_xticks(range(1,len(mons)+1))
                ax.set_xlim(0.5,len(mons)+0.5)
                ax.set_xticklabels(mons,rotation=90.)
                ax.text(0.02,0.98,reg,ha='left',va='top',transform=ax.transAxes)
                for run,col,mark in zip(runKeys,summaryColors,summaryMarkers):
                    D = np.zeros([len(mons),])
                    for a,mon in enumerate(mons):
                        if mon not in TS[reg][run].keys():
                            D[a] = np.nan
                            continue
                        ts = TS[reg][run][mon]
                        t = ts['model'][:,0]
                        N = float(len(t[np.isfinite(t)]))
                        for src in ['obs','model']:
                            ts[src][ts[src] > 100.] = np.nan

                        if norm == 'unit':
                            Tdist = ts['model'][:,0] - ts['obs'][:,0]
                            Sdist = ts['model'][:,1] - ts['obs'][:,1]
                        elif norm == 'max':
                            Tdist = (ts['model'][:,0] - ts['obs'][:,0]) / (np.nanmax(ts['obs'][:,0]) - np.nanmin(ts['obs'][:,0]))
                            Sdist = ts['model'][:,1] - ts['obs'][:,1] / (np.nanmax(ts['obs'][:,1]) - np.nanmin(ts['obs'][:,1]))
                        D[a] = np.nansum(np.sqrt(Tdist**2. + Sdist**2.)) / N

                    ax.scatter(np.arange(1,len(mons)+1),D,color=col,marker=mark)
                    maxD = np.max(np.array([maxD,np.nanmax(D)]))
                    minD = np.min(np.array([minD,np.nanmin(D)]))


            labels = []
            handles = []
            for run,col,mark in zip(runKeys,summaryColors,summaryMarkers):
                labels.append(run)
                handles.append(Line2D([],[],ls='none',color=col,marker=mark))
            fig.legend(handles,labels,bbox_to_anchor=(1.13,0.5))

            fac = 0.1
            for ax in axes[:,0]:
                ax.set_ylabel(u'$D_{TS}$/N')
                ax.set_ylim((1. - fac) * minD,(1. + fac) * maxD)
                if norm == 'max':
                    ax.set_yscale('log')
                    ax.set_ylim(10.**np.log10(minD),10.**np.log10(maxD))

            plt.tight_layout()

            os.makedirs(saveFolder,exist_ok=True)
            savePath = os.path.join(saveFolder,'Distance_TS_regional_summary.png')
            plt.savefig(savePath,dpi=300,bbox_inches='tight')
            plt.close()

        def plot_regional_TS_distance_individual(self,saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit'):
            mons = utils().get_monrange_TS(TS)
            regs = [k for k in sorted(coords.keys())]

            for reg in regs:
                print('Plotting TS distances for %s'%reg)
                fig,ax = plt.subplots(figsize=(7,5))
                ax.set_xticks(range(1,len(mons)+1))
                ax.set_xlim(0.5,len(mons)+0.5)
                ax.set_xticklabels(mons,rotation=90.)
                ax.grid(which='both')
                ax.text(0.02,0.98,reg,ha='left',va='top',transform=ax.transAxes)
                for run,col,mark in zip(runKeys,summaryColors,summaryMarkers):
                    D = np.zeros([len(mons),])
                    for a,mon in enumerate(mons):
                        if mon not in TS[reg][run].keys():
                            D[a] = np.nan
                            continue
                        ts = TS[reg][run][mon]
                        t = ts['model'][:,0]
                        N = float(len(t[np.isfinite(t)]))
                        for src in ['obs','model']:
                            ts[src][ts[src] > 100.] = np.nan

                        if norm == 'unit':
                            Tdist = ts['model'][:,0] - ts['obs'][:,0]
                            Sdist = ts['model'][:,1] - ts['obs'][:,1]
                        elif norm == 'max':
                            Tdist = (ts['model'][:,0] - ts['obs'][:,0]) / (np.nanmax(ts['obs'][:,0]) - np.nanmin(ts['obs'][:,0]))
                            Sdist = ts['model'][:,1] - ts['obs'][:,1] / (np.nanmax(ts['obs'][:,1]) - np.nanmin(ts['obs'][:,1]))
                        D[a] = np.nansum(np.sqrt(Tdist**2. + Sdist**2.)) / N

                    ax.scatter(np.arange(1,len(mons)+1),D,color=col,marker=mark)

                labels = []
                handles = []
                for run,col,mark in zip(runKeys,summaryColors,summaryMarkers):
                    labels.append(run)
                    handles.append(Line2D([],[],ls='none',color=col,marker=mark))
                ax.legend(handles,labels)

                ax.set_ylabel(u'$D_{TS}$/N')

                plt.tight_layout()

                os.makedirs(saveFolder,exist_ok=True)
                svPath = os.path.join(saveFolder,'Distance_TS_regional_%s.png'%reg)
                plt.savefig(svPath,dpi=300,bbox_inches='tight')
                plt.close()

        def plot_total_TS_distance(self,saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit'):
            mons = utils().get_monrange_TS(TS)
            regs = [k for k in sorted(coords.keys())]

            fig,ax = plt.subplots(figsize=(7,5))
            ax.set_xticks(range(1,len(mons)+1))
            ax.set_xticklabels(mons,rotation=90.)
            ax.set_xlim(0.5,len(mons)+0.5)
            ax.set_ylabel(u'$D_{TS}$/N')
            ax.grid(which='both')

            labels = []
            handles = []
            for run,col,mark in zip(runKeys,summaryColors,summaryMarkers):
                labels.append(run)
                handles.append(Line2D([],[],ls='none',color=col,marker=mark))

                D = np.zeros([len(mons),])
                N = np.zeros([len(mons),])
                for reg in regs:
                    for a,mon in enumerate(mons):
                        if mon not in TS[reg][run].keys():
                            continue
                        ts = TS[reg][run][mon]
                        t = ts['model'][:,0]
                        N[a] = N[a] + float(len(t[np.isfinite(t)]))
                        for src in ['obs','model']:
                            ts[src][ts[src] > 100.] = np.nan

                        if norm == 'unit':
                            Tdist = ts['model'][:,0] - ts['obs'][:,0]
                            Sdist = ts['model'][:,1] - ts['obs'][:,1]
                        elif norm == 'max':
                            Tdist = (ts['model'][:,0] - ts['obs'][:,0]) / (np.nanmax(ts['obs'][:,0]) - np.nanmin(ts['obs'][:,0]))
                            Sdist = ts['model'][:,1] - ts['obs'][:,1] / (np.nanmax(ts['obs'][:,1]) - np.nanmin(ts['obs'][:,1]))

                        D[a] = D[a] + np.nansum(np.sqrt(Tdist**2. + Sdist**2.))
                ax.scatter(np.arange(1,len(mons)+1),D/N,color=col,marker=mark)

            ax.legend(handles,labels)
            plt.tight_layout()

            os.makedirs(saveFolder,exist_ok=True)
            svPath = os.path.join(saveFolder,'Distance_TS_absolute.png')
            plt.savefig(svPath,dpi=300,bbox_inches='tight')
            plt.close()

    def plot_TS_distance(self,paths,metrics,coords,runKeys,summaryColors,summaryMarkers):
        print('Sorting cast data...')
        TS = utils().make_TS_dict(metrics,coords,runKeys)

        basePath = os.path.join(paths['figure_output_base_folder'],'TS_Diagrams','Monthly')
        print('Plotting monthly TS diagrams')
        plotters().TS_distance().plot_monthly_TS_diagrams(basePath,TS,coords,runKeys,summaryColors)
        print(' ')

        basePath = os.path.join(paths['figure_output_base_folder'],'TS_Diagrams','Annual')
        print('Plotting annual TS diagrams. Saving to %s'%basePath)
        plotters().TS_distance().plot_annual_TS_diagrams(basePath,TS,coords,runKeys,summaryColors)
        print(' ')

        saveFolder = os.path.join(paths['figure_output_base_folder'],'TS_Diagrams','TS_Distances_regional')
        print('Plotting regional subsets of TS distance. Saving to %s'%saveFolder)
        plotters().TS_distance().plot_regional_TS_distance_individual(saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit')
        print(' ')

        saveFolder = os.path.join(paths['figure_output_base_folder'],'TS_Diagrams')
        print('Plotting summary of regional subsets of TS distance. Saving to %s'%saveFolder)
        plotters().TS_distance().plot_regional_TS_distance(saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit')

        print('Plotting summary of total TS distance. Saving to %s'%saveFolder)
        plotters().TS_distance().plot_total_TS_distance(saveFolder,TS,coords,runKeys,summaryColors,summaryMarkers,norm='unit')
