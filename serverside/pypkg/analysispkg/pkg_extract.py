""" Helper functions for extractions
"""
import datetime as dt
import numpy as np

from analysispkg import pkg_data


def load_time(ncf, tstart, tend, time_offset=0, match='time', time_var=None):
    """
    Helper function to load time vector for variable varname and truncate to
    the range [tstart, tend].

    Parameters
    ----------
    ncf : netCDF4.Dataset object
        File to read time from.
    tstart,tend : datetime
        Time range to trim the extracted time series to.
    time_offset : integer
        Number of seconds to add to the time vector
    match : {'time','season'}, optional
        Option to match times:
            'time' matches exact times, should be used when comparing time series;
            'season' ignores year when matching times, should be used for constituents comparison.
        Default is 'time'.
    time_var : str, optional
        Name of the time variable in the file. Required to read time from model output files.
    Returns
    -------
        tuple: (truncation indices, truncated time vector)
    """
    timem1 = pkg_data.read_time(ncf, timename=time_var, offset=time_offset)
    if match == 'time':
        c, timem = match_time(timem1, tstart, tend)
    elif match == 'season':
        # c, timem = match_season(timem1, tstart, tend)
        #HACK To extract entire model series for HTG
        #TODO Didn't we want to extract entire series for time series variable in any case?
        c = np.arange(len(timem1))
        timem = timem1
    else:
        raise RuntimeError('Unknown option.')
    return c, timem

def match_time(timem1, tstart, tend):
    """ Helper for load_time()
    """
    c = np.where(np.logical_and(timem1 >= tstart, timem1 <= tend))[0]
    if len(c) == 0:
        return None, None
    else:
        return c, timem1[c]

def yearday(time):
    """ Decimal year-day.
    """
    timeyd = [t.timetuple().tm_yday +
              (t - t.replace(hour=0, minute=0, second=0)).total_seconds() / 3600 / 24
              for t in time]
    return np.array(timeyd)


def get_period_bounds(opt):
    """
    Simple helper function to find the limits of the analysis periods defined in opt
    """
    dates = []
    for period, (ana_start, ana_end) in opt['analysis']['periods'].items():
        dates += [ana_start, ana_end]
    return min(dates), max(dates)