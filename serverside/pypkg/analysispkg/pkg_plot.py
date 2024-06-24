# Various plotting utilities
import cmocean
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as clrs
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import scipy.stats as stats
from matplotlib.lines import Line2D

from analysispkg import pkg_data
from analysispkg import pkg_geo
from analysispkg import pkg_grid
from analysispkg import pkg_utils

# Go additions for bathy plot w/ rotation
from mpl_toolkits.basemap import Basemap
from matplotlib import patches


def plot_ticks(ax, minval=0, maxval=1, interval=0.1):
    """Plot utility: set the axis limits and tick interval according to the limits given by min/maxval"""

    upper = pkg_utils.fceil(maxval, roundto=interval)
    lower = pkg_utils.ffloor(minval, roundto=interval)

    if ax == 'x':
        plt.xlim(lower, upper)
        plt.xticks(np.arange(lower, upper, step=interval))
    elif ax == 'y':
        plt.ylim(lower, upper)
        plt.yticks(np.arange(lower, upper, step=interval))

    else:
        raise Exception('Unrecognized axis designation ' + ax)

    # ax.set_lim(lower,upper)
    # ax.set_ticks(np.arange(lower,upper,step=interval)


def axes_lims_with_margins(glam, gphi, data_lon, data_lat):
    data_lon = np.atleast_1d(data_lon)
    data_lat = np.atleast_1d(data_lat)
    # metrics and margins
    # diagonal extent of a grid cell
    dkm = pkg_data.haversine(glam[0, 0], gphi[0, 0], glam[1, 1], gphi[1, 1])
    ncells = 10  # approximate margin width around the data (in grid cells)
    cos_lat = np.cos(np.radians(data_lat[0]))
    lon_margin = ncells * dkm / (cos_lat * 111.325)
    lat_margin = ncells * dkm / 111.325
    map_xlim = [data_lon.min() - lon_margin, data_lon.max() + lon_margin]
    map_ylim = [data_lat.min() - lat_margin, data_lat.max() + lat_margin]
    return map_xlim, map_ylim, cos_lat


def plot_bathy_with_zoom(axm, bathy, data):
    """ Plot bathymetry for entire grid and an inset zoomed at the data

    Parameters
    ----------
    axm : list
        Length 2 list of axes for the main map and the inset
    bathy : dict
        Keys: lon, lat, depth
    data : dict
        Keys: lon, lat
    """
    # zoomed map limits
    map_xlim, map_ylim, _ = axes_lims_with_margins(bathy['lon'], bathy['lat'], data['lon'], data['lat'])

    # entire grid
    m1 = axm[0]
    _ = plot_bathy(m1, bathy, xlim=None, ylim=None, label_axes=False)
    # show zoom box in the entire grid plot
    m1.plot(np.array(map_xlim)[[0, 1, 1, 0, 0]], np.array(map_ylim)[[0, 0, 1, 1, 0]], 'k')

    # zoom on data area
    m2 = axm[1]
    b2 = plot_bathy(m2, bathy, xlim=map_xlim, ylim=map_ylim, label_axes=False)
    m2.plot(data['lon'], data['lat'], '.k', markersize=8)
    # cbar = plt.colorbar(b2, ax=m2, orientation="horizontal", pad=0.1)
    # cbar.ax.set_xlabel('Depth (m)')
    cbar = plt.colorbar(b2, ax=m2)
    cbar.ax.set_ylabel('Depth (m)')

    # disable scientific notation and offset for the lat,lon labels, and rotate the xlabels to avoid overlapping
    for ax in axm:
        ax.ticklabel_format(useOffset=False, style='plain')
        for tick in ax.get_xticklabels():
            tick.set_rotation(20)
            #tick.set_ha('right')


def get_bathy(opt):
    lon, lat = pkg_grid.tgrid(opt['file_coord'])
    lonf, latf = pkg_grid.fgrid(opt['file_coord'])
    lonfe, latfe = pkg_geo.expandf(lonf,latf) # corner grid
    depth = pkg_grid.bathymetry(opt['file_bathy'], maskfile=opt['file_mesh'])
    bathy = dict(lon=lon, lat=lat, depth=depth, lonfe=lonfe, latfe=latfe)
    return bathy


def plot_bathy(ax, bathy, land_color='grey', xlim=None, ylim=None, label_axes=False, vmin=None, vmax=None):
    """ Plots bathymetry as pcolormesh in axes `ax` """
    cmap = cmocean.cm.deep
    cmap.set_bad(color=land_color)
    h = ax.pcolormesh(bathy['lonfe'], bathy['latfe'], bathy['depth'], cmap=cmap, vmin=vmin, vmax=vmax,
                      shading='auto', rasterized=True)
    if xlim is not None and ylim is not None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        mean_lat = (ylim[0] + ylim[1]) / 2
    else:
        mean_lat = (bathy['lat'].min() + bathy['lat'].max()) / 2
    if label_axes:
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    ax.set_aspect(1 / np.cos(np.radians(mean_lat)))
    return h



########################################################################

# try to plot with rotation -GO change 202309
#def plot_bathy_GO(ax, bathy, 
#                  land_color='grey', 
#                  xlim=None, ylim=None, 
#                  label_axes=False, 
#                  vmin=None, 
#                  vmax=None):
#                  
#    """ Plots bathymetry as pcolormesh in axes `ax` """
#    cmap = cmocean.cm.deep
#    cmap.set_bad(color=land_color)
#    h = ax.pcolormesh(bathy['lonfe'], bathy['latfe'], bathy['depth'], cmap=cmap, vmin=vmin, vmax=vmax,
#                      shading='auto', rasterized=True)
#    if xlim is not None and ylim is not None:
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#        mean_lat = (ylim[0] + ylim[1]) / 2
#    else:
#        mean_lat = (bathy['lat'].min() + bathy['lat'].max()) / 2
#    if label_axes:
#        ax.set_xlabel('Longitude')
#        ax.set_ylabel('Latitude')
#    ax.set_aspect(1 / np.cos(np.radians(mean_lat)))
#    return h

def adjust_map_GO(ax, m,
               lat_bl, 
               lon_bl, 
               lat_br, 
               lon_br,
               lat_tl,
               lon_tl,
               lat_bl2,
               lon_bl2
              ):
    
    # set width using map units
    # bottom left
    x_bl,y_bl = m(lon_bl,lat_bl)
    x_br,_ = m(lon_br,lat_br)
    ax.set_xlim(x_bl, x_br)

    # top left
    x_tl, y_tl = m(lon_tl, lat_tl)
    x_bl,y_bl = m(lon_bl2,lat_bl2)
    ax.set_ylim(y_bl, y_tl)

    # fix a little path in bottom right
    lccx_TL, lccy_TL = m(-122.83, 49.4)
    lccx_BR, lccy_BR = m(-122.58, 48.7)
    lccx_BL, lccy_BL = m(-122.33, 48.7)
    lccw = lccx_BL - lccx_BR
    lcch = lccy_TL - lccy_BL

    ax.add_patch(patches.Rectangle(
                (lccx_BL, lccy_BL), lccw, lcch, 
                facecolor='burlywood', edgecolor='k',
                linewidth=0,
                zorder=0))


def plot_bathy_GO(ax, grid, w_map=[-124, -123.9, 47.7, 50.6], 
             rotation=39.2, 
             par_inc=0.25,
             mer_inc=0.5, 
             fs=7,
             bg_color='#969696'
            ):
    """
    """

    # Make projection
    m = Basemap(ax=ax, 
                projection='lcc', resolution='c',
                lon_0=(w_map[1] - w_map[0]) / 2 + w_map[0] + rotation,
                lat_0=(w_map[3] - w_map[2]) / 2 + w_map[2],
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3])
    
    land_color = 'burlywood'
    cmap = cmocean.cm.deep
    cmap.set_bad(color=land_color)
    x, y = m(grid['lonfe'], grid['latfe'])
    h = ax.pcolormesh(x, y, grid['depth'], cmap=cmap, #vmin=vmin, vmax=vmax,
                      shading='auto', rasterized=True)

    m.drawmeridians(np.arange(-125.5, -122, mer_inc), labels=[0, 0, 0, 0], linewidth=0.2, fontsize=fs)
#     m.drawparallels(np.arange(48, 51, par_inc), labels=[1, 0, 0, 0], linewidth=0.2, fontsize=fs)
    m.drawparallels(np.arange(47, 51, par_inc), labels=[0, 0, 0, 0], linewidth=0.2, fontsize=fs)   
    
    # Add features and labels
    x, y = m(grid['lon'], grid['lat'])
    
    ax.contourf(x, y, grid['depth'], [-0.01, 0.01], colors=land_color)
    ax.contour(x, y, grid['depth'], [-0.01, 0.01], colors='black', linewidths=0.1)
    #ax.contourf(x, y, grid['depth'], [0.011,500], colors=bg_color)
#     m.drawmeridians(np.arange(-125.5, -122, mer_inc), labels=[0, 0, 0, 1], linewidth=0.2, fontsize=fs)


    return m
    

########################################################################  
    
    




def plot_map(ax, lm, xlim, ylim, clr, label_axes=False):
    """Takes a landmask object and plots it in axes ax, with lat/lon limits set by lim."""

    c = clrs.ListedColormap([clr, clr])
    ax.pcolormesh(lm.lon, lm.lat, lm.land, edgecolors=None, cmap=c, shading='auto')  # may be better as contourf?
    # ax.pcolormesh(lm.lon, lm.lat, lm.land, edgecolors=None, cmap=c, shading='auto', rasterized=True)  # may be better as contourf?
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if label_axes:
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')


# GO changed args to include linewidth
def plot_lcn(ax, loc, lbl, clr, sz=60, lw=0.6):
    """Takes a tuple with the latlon coordinates of a point, and plots a dot of size sz at that point."""
    #ax.scatter(loc[0], loc[1], c=clr, edgecolors='none', s=sz, label=lbl)
    # GO change to change location to empty circle
    ax.scatter(loc[0], loc[1], c='none', edgecolors=clr, s=sz, linewidth=lw, label=lbl)


# added by GO 20231003 - for plotting a set of stations as numbers on map
def plot_lcn_GO(ax, loc, lbl, clr, sz=40):
    """
    Takes a list of tuples with latlon coordinates (locs) and a list of labels (lbls),
    and plots dots with corresponding number symbols on the given axis (ax).
    """
    if len(lbl) > 1:
      sz = sz * 2
    ax.scatter(loc[0], loc[1], c=clr, edgecolors='none', s=sz, marker=f'${lbl}$')



def plot_lcn_transect(ax, locs, lbl, clr, lw=2):
    """Takes a time series of locations, draws lines between them."""

    ax.plot(locs[0], locs[1], color=clr, lw=lw, label=lbl)


def list_info(info, dictionary, title="", config=""):
    """ Uses an otherwise blank pair of axes to reserve space, and prints metadata about the station / data.
    Info is an indeterminate list of tuple pairs.  First element is the key in the dictionary to the string denoting
    the type of info being printed. """

    txtstr = ""

    if len(title) > 0:
        txtstr = title + '\n'
    if len(config) > 0:
        txtstr = txtstr + config + '\n'

    for i in info:
        h = dictionary[i[0]]
        if len(h) > 0:
            txtstr = txtstr + h + ": " + str(i[1]) + '\n'
        else:
            txtstr = txtstr + str(i[1]) + '\n'

    return txtstr


def axesetc(ax, pltopt, timeaxis, xlabel, ylabel):
    """Subroutine for the bog standard axes plotting: add a visible zero, grid lines, axes labels, etc."""

    # TODO: Unhardwire this?
    # ax.plot(timeaxis, np.zeros(timeaxis.shape), color=pltopt['colors']['zero'], linewidth=pltopt['residual']['lw'])
    # ax.plot(timeaxis, np.zeros(timeaxis.shape), color=pltopt['colors']['zero'], linewidth=2)
    # TODO Refactor to exclude `timeaxis` argument
    ax.axhline(0., color=pltopt['colors']['zero'], linewidth=2)  # , ls='--'
    ax.grid(color=pltopt['colors']['grid'], linestyle=':', linewidth=1)
    # ax.grid(color='#888888', linestyle=':', linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)  # a1.legend()


def set_colourbar_evenly():
    """ For difference plots, ensure that there's an equal amount of colourbar on either side so visual comparison
    is fair. """
    vmin, vmax = plt.gci().get_clim()
    vamp = np.nanmax(np.absolute([vmin, vmax]))
    plt.gci().set_clim(-vamp, vamp)


# Subroutine from https://gist.github.com/jakevdp/91077b0cae40f8f8244a
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
    # return LinearSegmentedColormap.from_list(cmap_name, color_list, N)


# set the colormap and centre the colorbar
# http://chris35wills.github.io/matplotlib_diverging_colorbar/
class MidpointNormalize(clrs.Normalize):
    """
        Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

        e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
        """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        clrs.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def get_eval_str(where, ismap=False):
    """ Returns a grammatically sensible string to use in plot title indicating whether plot show x metres above
    the seabed or below the surface. """

    if ismap is False:
        if where == 'bottom':
            refstr = 'above sea bed'
        elif where == 'surface':
            refstr = 'below surface'
        else:
            refstr = None
    else:
        if where == 'bottom':
            refstr = 'Closest to Sea Bed'
        elif where == 'surface':
            refstr = 'Nearest the Surface'
        else:
            refstr = None
    return refstr


def figure_saver(filename_without_suffix, file_formats=None, fig=None, **options):
    """ Helper for figure saving.
    It should be used for all figure saving instead of plt.savefig
    Parameters
    ----------
    filename_without_suffix : str
    file_formats : list
    fig : <matplotlib.figure>, optional
        Default is current figure
    **options : kwargs
        Options for plt.savefig()
    """
    dname, fname = os.path.split(filename_without_suffix)
    os.makedirs(dname, exist_ok=True)
    filename_without_suffix_clean = os.path.join(dname, fname.replace("'","").replace(" ","_"))
    if file_formats is None:
        file_formats = ['png']
    if fig is None:
        fig = plt.gcf()
    for file_format in file_formats:
        fig.savefig(f'{filename_without_suffix_clean}.{file_format}', format=file_format, **options)
    if len(file_formats) > 1:
        msg = "Saved figure {}.{{{}}}".format(filename_without_suffix_clean,','.join(file_formats))
    else:
        msg = "Saved figure {}.{}".format(filename_without_suffix_clean, file_formats[0])
    return msg


# Common SST, LH routines
def get_analysis_data(opt, analysis_files, casenames=None):
    data = {} 
    for f in analysis_files:
        with open(f,'rb') as fid:
            newdat = pickle.load(fid)
        for k in newdat.keys():
            if k not in data.keys():
                data[k] = newdat[k]
            # Need the full gamut of observations available, even if not available in the first model loaded
            if k == 'obs':
                data[k].update(newdat[k])

    # For backward compatability with single plotting routines, which don't have a casenames string quickly available
    if casenames is not None:    
        stationlists = []
        for c in casenames:
            if c in data.keys():
                stationlists.append(list(data[c].keys()))
            else:
                stationlists.append([])

        # Determine which stations are available in all sources (obs, reference model and candidate models)
        # Now with union!
        stations = pkg_utils.get_common_stations(opt, casenames, stationlists)

    else:   # Keep the original block, which is probably not necessary for single config plotting
        stations = None
        inst = None
        for src in data.keys():
            if inst is None:
                inst = set(data[src].keys())
            else:
                inst = inst.intersection(set(data[src].keys()))

        # Remove stations not in this common list
        for src in data.keys():
            keys = list(data[src].keys())
            for s in keys:
                if s not in inst:
                    data[src].pop(s)

    return data, stations


def get_plot_params(opt, mods=None):
    if mods is None:
        # We're doing single-case plots
        reference = opt['casename']
        candidates = []
        colours = {'obs': opt['plot']['colors']['obs'],
                   'reference': opt['plot']['colors']['mod'],
                   'candidates': ['']}
        markers = {'obs': 'x',
                   'reference': 'o',
                   'candidates': ['']}
        lineweights = {'obs': 2, 'reference': 1, 'candidates': [1]}
        formats = opt['plot']['file_formats']

    else:
        # We're doing intercomparison plots
        reference = mods[0]['casename']
        candidates = [x['casename'] for x in mods[1:]]

        colours = {'obs': 'k',
                   'reference': opt['mods'][0]['color'],
                   'candidates': [x['color'] for x in opt['mods'][1:]]}
        markers = {'obs': 'x',
                   'reference': opt['mods'][0]['marker'],
                   'candidates': [x['marker'] for x in opt['mods'][1:]]}
        lineweights = {'obs': 2,
                       'reference': opt['mods'][0]['lineweight'],
                       'candidates': [x['lineweight'] for x in opt['mods'][1:]]}
        formats = opt['file_formats']
    return reference, candidates, colours, markers, lineweights, formats


def populate_statistics(statistics, v, lab, rows):
    if np.isnan(v).all():
        return statistics, rows

    statistics['mean'][lab] = np.nanmean(v)
    # GO edit  May 26 2023
    #statistics['rms'][lab] = np.sqrt(np.nanmean(v) ** 2.)
    statistics['rms'][lab] = np.nanstd(v)
    statistics['skew'][lab] = stats.skew(v[~np.isnan(v)])
    statistics['kurt'][lab] = stats.kurtosis(v[~np.isnan(v)], fisher=False)

    rows.append(lab)

    return statistics, rows


def plot_instrument(data, s, reference, candidates, colors, lineweights, tsname, svFolder='./', file_formats=None,
                    display=False):
    statistics = {}; rows = []
    for stat in ['mean','rms','skew','kurt']:
        statistics[stat] = {}

    #TODO Add to parameters
    # start = datetime.datetime(2017,7,21)
    # end = datetime.datetime(2018,1,1)
#        end = datetime.datetime(2017,8,3)
    start = data['obs'][s]['start']
    end   = data['obs'][s]['end']

    data_obs = data['obs'][s]
    timeF = data_obs['filt_interp_data']['time']
    c = np.where(np.logical_and(timeF >= start,timeF <= end))[0]
    time = timeF[c]
    To = data_obs['filt_interp_data'][tsname][c]

    statistics, rows = populate_statistics(statistics, To, 'Observation', rows)

    fig = plt.figure(figsize=(10.5,8.))
    gs = gridspec.GridSpec(3,4)

    axDT = fig.add_subplot(gs[1,:-1])
    axDTh = fig.add_subplot(gs[1,-1])
    axT = fig.add_subplot(gs[0,:-1],sharex=axDT)
    axTh = fig.add_subplot(gs[0,-1])

    axTextT = fig.add_subplot(gs[2,:])

    for ax in [axDT,axDTh]:
        ax.axhline(0.,color='k',ls='--')

    plt.sca(axT)
    plt.tick_params(axis='x',which='both',labelbottom=False)

    plt.sca(axTh)
    plt.tick_params(axis='x',which='both',bottom=False,top=True,labelbottom=False,labeltop=True)
    axTh.xaxis.set_label_position("top")
    axTh.xaxis.tick_top()

    for ax in [axTh,axDTh]:
        plt.sca(ax)
        plt.tick_params(axis='y',which='both',left=False,right=True,labelleft=False,labelright=True)

    for ax in [axTh,axDTh]:
        ax.set_xlabel('Probability Density')

    axT.plot(time,To,color=colors['obs'],alpha=0.7,lw=lineweights['obs'],zorder=1001)

    if not np.all(np.isnan(To)):
        axTh.hist(To,color=colors['obs'],histtype='step',orientation='horizontal',density=True,lw=lineweights['obs'],zorder=1001)

    for ax in [axDT,axDTh,axT,axTh]:
        ax.grid(which='both')

    for ax in [axDT,axT]:
        ax.set_xlim(start,end)

    if tsname == 'temperature':
        axT.set_ylabel('SST (degC)')
        axDT.set_ylabel(r'$\Delta$T (degC)' + '\nModel - Observations')
        abbr, deltachar = 'SST','T'
    elif tsname == 'salinity':
        axT.set_ylabel('SSS (PSU)')
        axDT.set_ylabel(r'$\Delta$S (PSU)' + '\nModel - Observations')
        abbr, deltachar = 'SSS','S'

    for mod,col,LW in zip([reference] + candidates,\
                          [colors['reference']] + colors['candidates'][:len(candidates)],\
                          [lineweights['reference']] + lineweights['candidates'][:len(candidates)]):

        timeM = data[mod][s]['filt_interp_data']['time']
        c = np.where(np.logical_and(timeM >= start, timeM <= end))[0]
        timeM = timeM[c]
        Tm = data[mod][s]['filt_interp_data'][tsname][c]
        if np.isnan(Tm).all():
            continue

        # Model and obs may not be on the same time vector!
        # This is because we analyze models separately and the time interpolation in SST analysis code was
        # designed to analyse all of the models in a comparison at once. We need to rethink this somehow to take
        # into account that models may have different output frequencies (daily vs hourly)...
        # -- Interim workaround for now:  interpolate onto obs times here
        tm = np.array([(x - start).total_seconds() for x in timeM])
        to = np.array([(x - start).total_seconds() for x in time])
        Tm = np.interp(to, tm, Tm)
        # -- End interim workaround

        statistics, rows = populate_statistics(statistics, Tm, mod, rows)

        dT = Tm - To

        statistics, rows = populate_statistics(statistics, dT,
                                               r'$\Delta {%s}_{%s}$' % (deltachar, mod.replace('_', '-')), rows)

        axDT.plot(time,dT,color=col,lw=LW,zorder=1000)
        axDT.tick_params(axis='x', labelrotation=45)

        if not np.all(np.isnan(dT)):
            axDTh.hist(dT,color=col,lw=LW,histtype='step',orientation='horizontal',density=True,zorder=1000)

        axT.plot(time,Tm,color=col,lw=LW,zorder=1000)

        if not np.all(np.isnan(Tm)):
            axTh.hist(Tm,color=col,lw=LW,histtype='step',orientation='horizontal',density=True,zorder=1000)

    columns = ['Mean/Bias','St.Dev./RMS','Skewness','Kurtosis']

    table_data = np.zeros([len(rows),len(statistics)]).astype(str)
    for j,src in enumerate(rows):
        for i,stat in enumerate(statistics.keys()):
            table_data[j,i] = str(np.round(statistics[stat][src],3))

    axTextT.table(cellText=table_data,bbox=[0.2,0.,0.8,1.],colLabels=columns,rowLabels=rows,fontsize=8)
    axTextT.axis('off')
    axTextT.set_title(abbr + ' Statistics')

    handles = [Line2D([],[],color=col,lw=LW) for col,LW in zip([colors['obs']] + [colors['reference']] + \
                                                               colors['candidates'][:len(candidates)],\
                                                               [lineweights['obs']] + [lineweights['reference']] + \
                                                               lineweights['candidates'][:len(candidates)])]
    labels = ['Observations'] + [reference] + candidates
    #axT.legend(handles, labels, bbox_to_anchor=(1., 0.2))
    axT.legend(handles,labels,ncol=len(labels),bbox_to_anchor=(1.0,-0.1))
    axT.set_title("Station " + s) # GO changed - LH not 'buoys'

    plt.tight_layout()

    if not display:
        os.makedirs(svFolder,exist_ok=True)
        filename_without_suffix = os.path.join(svFolder, f'{abbr}_{s}')
        figure_saver(filename_without_suffix, file_formats, fig, dpi=450, bbox_inches='tight')
        plt.close()
