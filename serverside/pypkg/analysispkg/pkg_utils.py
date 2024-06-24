import argparse
import collections.abc
import datetime as dt
import functools
import getpass
import json
import logging
import math
import numpy as np
import os
import pdb
import sys
import time
import yaml

INSTRUMENTS = ['CTD', 'LH', 'SST', 'OBS']


debug_switch=[False]
def debug(f):
    """
    Decorator to launch pdb on an unhandled exception if debug_switch is True
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        if not debug_switch[0]:
            return f(*args, **kwargs)
        try:
            print("pkg_utils.debug activated for",f)
            return f(*args, **kwargs)
        except Exception:
            pdb.post_mortem(sys.exc_info()[2])
    return wrapper


def positive_int(inval):
    val = int(inval)
    if val >= 1:
        return val
    else:
        raise argparse.ArgumentTypeError('Value must be >= 1.')

def try_and_log(func):
    """ Decorator to log the error in case a function fails.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            out = func(*args, **kwargs)
        except Exception as e:
            out = None
            logging.error(e, exc_info=True)
        return out
    return wrapper


def time_it(func):
    """ Decorator to time a function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        out = func(*args, **kwargs)
        print("--- {} run time: {:.2f} seconds ---" .format(func.__name__,
              round(time.time() - start_time, 2)))
        return out
    return wrapper


def parse_time (str_date):
    """Utility: Converts an eight digit data into year, month, and day strings"""

    if str_date.isdigit:

        str_yr = str_date[0:4]
        str_mth = str_date[4:6]
        str_day = str_date[6:8]
    else:
        raise Exception('{}  contains non-digit characters.'.format(str_date))


    return str_yr, str_mth, str_day


def numdate(str_date):
    """Utility: Takes a formatted date and returns a datetime date."""
    fmt_ymd = '%Y%m%d'
    fmt_ymd_hms = '%Y%m%d %H:%M:%s'
    if len(str_date)==8:
        num_date = dt.datetime.strptime(str_date, fmt_ymd)
        #num_date = calendar.timegm(time.strptime(str_date, '%Y%m%d'))
    elif len(str_date)==17:
        num_date = dt.datetime.strptime(str_date, fmt_ymd_hms)
    else:
        raise Exception('Date format must be either {} or {}'.format(fmt_ymd,fmt_ymd_hms))
    return num_date

def strdate(dt_date):
    """Utility: Takes a datetime date and returns a string."""
    fmt_ymd = '%Y%m%d'
    if isinstance(dt_date, dt.datetime):
        str_date = dt_date.strftime(fmt_ymd)
    else:
        raise Exception('strdate requires a datetime variable, not {}'.format(type(dt_date)) )
    return str_date



def time_range (startdate, enddate, interval, offset=0.0):
    """Utility: Generates a consistently spaced time array from startdate to enddate
                offset is used if required timestamps do not fall on the hour"""

    # dates are strings, and interval is in hours
    # convert to datetime object
    sd = numdate(startdate)
    ed = numdate(enddate)

    # get the number of total intervals (hourly, halfhourly, whatever) between the requested dates
    # use seconds so everything can be integers.
    window = int( (ed - sd).total_seconds() )
    nsec   = int (interval * 3600)

    # get an array of all the timestamps at the specified interval
    arr_time = np.asarray ([ sd + dt.timedelta(hours=float(offset)) + dt.timedelta(seconds=t) for t in range (0, window, nsec) ])

    return arr_time


# These three functions round a value to the nearest multiple of a float (can be < 1).  Combine?
def fceil (val, roundto=0.01, precision=2):	
    return round(roundto * math.ceil(float(val)/roundto),precision)
def ffloor (val, roundto=0.01, precision=2):	
    return round(roundto * math.floor(float(val)/roundto),precision)
def fround (val, roundto=0.01, precision=2):	
    return round(roundto * round(float(val)/roundto),precision)


def load_config_yaml(f):
    configpath = os.path.normpath(os.path.dirname(__file__) + "/../config")
    yamldata = load_yaml(os.path.join(configpath, f))
    return yamldata


def load_yaml(yamlfile):
    """ Helper to load a YAML
    """
    def date_to_datetime(loader, node):
        """ The default YAML loader interprets YYYY-MM-DD as a datetime.date object
            Here we override this with a datetime.datetime object with implicit h,m,s=0,0,0 """
        d = yaml.constructor.SafeConstructor.construct_yaml_timestamp(loader,node)
        if type(d) is dt.date:
            d = dt.datetime.combine(d, dt.time(0, 0, 0))
        return d
    yaml.constructor.SafeConstructor.yaml_constructors[u'tag:yaml.org,2002:timestamp'] = date_to_datetime
    with open(yamlfile, 'r') as ya:
        try:
            yamldata = yaml.safe_load(ya)
        except Exception as e:
            print("Error importing YAML file {} {})".format(yamlfile,e))
            raise
    return yamldata


def load_options_comparison(comparisonyaml):
    """ Loads default comparison options and override with a comparison options file.
    """
    # Import the reference options
    configpath = os.path.normpath(os.path.dirname(__file__) + "/../config")
    opts = load_yaml(os.path.join(configpath, "default_comparison.yaml"))

    # Import the comparison options
    opts_comparison = load_yaml(comparisonyaml)
    opts_comparison['optfile'] = comparisonyaml
    check_deprecations_comparison(comparisonyaml, opts_comparison)

    # Apply overrides
    opts = update(opts, opts_comparison)

    # Load model options for each and store in list
    mods = []
    for case in opts['mods']:
        mods += [ load_options(case['optfile']) ]

    # Decide which domain to use for the map icons
    # First try to take the first model that has "map: True",
    # and fallback to first model if no model with "map: True" found
    def getplotmod(opts):
        for i,cmod in enumerate(opts['mods']):
            if cmod.get('map',False):
                return i
        return 0
    p = getplotmod(opts)
    opts['plotmod'] = p

    # Populate opt with some variables from the first of the mods in the yaml.
    opts['file_coord'] = mods[p]['file_coord']
    opts['file_bathy'] = mods[p]['file_bathy']
    opts['file_mesh']  = mods[p]['file_mesh']

    opts['plot']['colors'] = mods[p]['plot']['colors']
    opts['plot']['fonts']  = mods[p]['plot']['fonts']

    opts['dir_plots'] = opts['output_path']
    opts['dir_logs'] = os.path.join(opts['output_path'],'LOGS')

    os.makedirs(opts['dir_plots'], exist_ok=True)

    return opts,mods


def load_options(yamlfile):
    """ Load default options and override with analysis options.
    Loads default options and updates from a model config (if specified) and
    the runid config (if specified), populates obs and sets directories.
    """
    # Import the reference options
    configpath = os.path.normpath(os.path.dirname(__file__) + "/../config")
    opts = load_yaml(os.path.join(configpath, "default.yaml"))

    # Import the analysis options
    opts_analysis = load_yaml(yamlfile)
    opts_analysis['optfile'] = yamlfile
    check_deprecations(yamlfile, opts_analysis)

    # Apply overrides
    opts = update(opts, opts_analysis)

    # Todo: consolidate these two into one
    opts['casename'] = opts['config'] + '-' + opts['runid']
    opts['configname'] = opts['config'] + '-' + opts['runid']

    # Generate directory structure
    opts = set_reference_paths(opts)

    # check paths
    checkpaths(opts, ['file_bathy', 'file_coord', 'file_mesh', 'src_output_dir', 'obs:root'])

    return opts

def checkpaths(opt,pathkeys):
    error=False
    for pathkey in pathkeys:
        val = opt
        for key in pathkey.split(':'):
            val = val[key]
        if val is not None and not os.path.exists(val):
            print("Input path {} {} does not exist, can not proceed".format(pathkey,val))
            error=True
    if error:
        raise FileNotFoundError

def check_deprecations(fname,opt):
    prefix = "check_deprecations for config file {}: ".format(fname)
    print('checking keys - GO')
    print(opt.keys()) # GO check 2024
    if 'TG' in opt.keys():
        raise ValueError(prefix + 'TG root should be split and moved to analysis:TG and plot:TG, see default.yaml for example')
    if 'plot_constituents' in opt.keys():
        raise ValueError(prefix + 'plot_constituents root is obsolete and should be moved to plot:TG:constituents')
    if 'list_constituents' in opt['plot'].keys():
        raise ValueError(prefix + 'map:list_constituents is obsolete and should be removed')
    if 'time_getter' in opt.keys():
        print(prefix + "WARNING: time_getter is no longer used and is ignored in your options file")
    if 'start_date' in opt.keys():
        print(prefix + "WARNING: start_date is no longer used and is ignored in your options file")
    if 'final_date' in opt.keys():
        print(prefix + "WARNING: final_date is no longer used and is ignored in your options file")

def check_deprecations_comparison(fname,opt):
    prefix = "PyAP comparison config file {}: ".format(fname)
    if 'TG' in opt.keys() and 'stations' in opt['TG'].keys():
        raise ValueError(prefix+'TG:stations is no longer used, stations are detected automatically')

def prologue(opt):
    print("PyAP starting up")
    print("Config file: {}".format(opt['optfile']))
    fname = printopts(opt)
    print("Saving combined (defaults + overrides) options to {}".format(fname))
    print("nproc: {}".format(opt['nproc']))
    print("----------")

def update(d, u):
    """ Recursive update for a nested dictionary. """
    # Source: https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        elif (type(d) is list and len(d) == 0) or d is None:
            d = {k: v}  # empty, no such entry in default.yaml
        else:
            d[k] = v
    return d


def printopts(opts):
    """Pretty prints the opts dict to a file.
       This is useful to confirm it set up the right paths and that overrides are working, etc.
    """
    os.makedirs(opts['output_path'], exist_ok=True)
    fname = os.path.join(opts['output_path'], 'runtime_opts.yaml')
    info1 = '# combined options file generated at {}\n'.format(str(dt.datetime.now()))
    info2 = '# command line: {}\n\n'.format(' '.join(sys.argv))
    opts_pp = yaml.dump(opts,default_flow_style=None)
    with open(fname,'w') as f:
        f.write(info1+info2)
        f.write(opts_pp)
    return fname


def getuservar(var):
    configpath = os.path.normpath(os.path.dirname(__file__) + "/../config")
    userpaths = load_yaml(os.path.join(configpath, "users.yaml"))
    username = getpass.getuser()
    try:
        value = userpaths[username][var]
    except:
        print("Error loading user {} variable {}, please check config/users.yaml".format(username,var))
        raise
    return value


def set_reference_paths(opt):
    """
    Sets derived paths according to the following structure:

    The output_dir is specified in the analysis yaml

    Extracted data: <runroot>/EXTRACT/INSTRUMENT
    Processed data: <runroot>/PROCESS/INSTRUMENT
    Plots etc:      <runroot>/PLOTS/INSTRUMENT
    """
    runroot = opt['output_path']
    opt['dir_run']     = runroot  # not sure if we need this.
    opt['dir_extract'] = os.path.join(runroot, 'EXTRACT')  # Model extractions
    opt['dir_process'] = os.path.join(runroot, 'PROCESS')  # Interpolated/filtered data, calculated scores, etc
    opt['dir_plots']   = os.path.join(runroot, 'PLOTS')    # Plots, tables, figures of all sorts
    opt['dir_report']  = os.path.join(runroot, 'REPORT')   # Auto-generated report
    opt['dir_logs']    = os.path.join(runroot, 'LOGS')     # Logs for submitted jobs

    return opt


def load_index(archive, instrument):
    """ Load the indexed data from the obs archive for 'instrument'
    """
    indexfile = os.path.join(archive['root'], 'INDEX', instrument + ".json")
    obsindex = []
    if os.path.exists(indexfile):
        with open(indexfile, "r") as fid:
            obsindex = json.load(fid,object_hook=__json_decode_date)
        # relative path to absolute path
        for inst in obsindex:
            inst['filename'] = os.path.join(archive['root'], inst['filename'])
    return obsindex


def __json_decode_date(dic):
    if 'time' in dic:
        if dic["time"] is str:
            dic["time"] = dt.datetime.fromisoformat(dic["time"])
        elif type(dic["time"]) is list:
            dic["time"] = [dt.datetime.fromisoformat(x) for x in dic["time"]]
        return dic


def get_process_dir(opt, instrument, period):
    """ Generate directory name and full path to the analysis file """
    if period is None:
        outpath = os.path.join(opt['dir_process'], instrument)
    else:
        outpath = os.path.join(opt['dir_process'], instrument, period)
    return outpath


def get_process_path(opt, instrument, period):
    """ Generate directory name and full path to the analysis file """
    outpath = get_process_dir(opt, instrument, period)
    filename = instrument + '_class4_' + opt['casename'] + ".pickle"
    filepath = os.path.join(outpath, filename)
    return outpath, filepath


def get_plot_dir(opt, instrument, period=None, where_eval=None):
    """ Generate directory name for plots/scores """
    outpath = os.path.join(opt['dir_plots'], instrument)
    if period is not None:
        outpath = os.path.join(outpath, period)
    return outpath


def csv_file(typ, period, scoretype, statistic, opt):
    """ Generate file name for scores """
    return '_'.join((typ, period, scoretype, statistic, opt['runid'] + '.csv'))


def parse_constituents (incons):
    if type(incons) is list:
        outcons = incons
    elif incons == 'cons_5':
        outcons = ['M2', 'K1', 'S2', 'N2', 'O1']
    else:
        outcons = None
    return outcons


def get_common_stations(opt, casenames, stationlists):
    """
    opt: comparisons options file
    casenames: list of casenames
    stationlists: list of list of stations, one list of stations per casename

    The baseline common stations is the set intersection across all models
    This set is expanded by set union with stations that have union:True specified in the comparison options
    """
    def union_enabled(opt, casename):
        """
        Helper to check if union:True is specified for this model
        """
        if opt.get('mods', None):  # For compatability with single-config setup, since this is used in ACOM routines that are common to both
            for mod in opt['mods']:
                modopt = load_options((mod['optfile']))
                if modopt['casename'] == casename:
                    return mod.get('union', False)
        return False
    # Reduce by set intersection
    r = set.intersection(*[set(x) for x in stationlists])
    # Expand by set union
    for stationlist, casename in zip(stationlists, casenames):
        if union_enabled(opt, casename):
            r = set.union(r, set(stationlist))
    return r
