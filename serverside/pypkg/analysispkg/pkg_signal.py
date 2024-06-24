""" Signal processing utilities
"""
import datetime
from collections import Counter

import numpy as np
import scipy.cluster
import scipy.linalg
import scipy.signal
import scipy.stats
import scipy.interpolate
from scipy.spatial import cKDTree

from analysispkg import pkg_interp


def interp(ti, t, v, left=None, right=None, period=None):
    # np.interp for times in datetime format
    # left, right : extrapolation values; extrapolates with nearest to both sides by default
    vi = np.interp(numtime(ti), numtime(t), v, left=left, right=right, period=period)
    return vi


def numtime(t):
    # Convert datetime to timestamp, float seconds since January 1, 1970 at UTC (on UNIX)
    # return np.array([x.timestamp() for x in t])
    return [x.timestamp() for x in t]


def numtime_to_datetime(tstamp):
    # Convert timestamp to datetime
    return np.array([datetime.datetime.fromtimestamp(x) for x in tstamp])


def regularize(t, v, si):
    """ Interpolate irregularly sampled series to regular intervals
    
    First time stamp in the output series coincides with the first time stamp
    in the original series.
    
    Parameters
    ----------
    
    t : (N,) array_like
        x-coordinate; can be array of datetime.
    v : (N,...) array_like
        corresponding y-values
    si : datetime.timedelta
        Target sampling interval, e.g. datetime.timedelta(minutes=2)

    Returns
    -------
    
    tr,vr : array-like
        Regularized series
    """
    tr = np.arange(t[0], t[-1] + si, si).astype(datetime.datetime)
    # vr = np.interp(matplotlib.dates.date2num(tr),matplotlib.dates.date2num(t),v)
    # remove times above the interpolation range
    tr = tr[tr <= t[-1]]
    f = scipy.interpolate.interp1d(numtime(t), v, axis=0, assume_sorted=True)
    vr = f(numtime(tr))
    return tr, vr


def regularize_num(tnum, v, interval):
    """ Interpolate irregularly sampled series to regular intervals

    First time stamp in the output series coincides with the first time stamp
    in the original series.

    Same as regularize(), but works with time as a float array.

    Parameters
    ----------

    tnum : (N,) array_like
        x-coordinate.
    v : (N,...) array_like
        corresponding y-values
    interval : float
        Target sampling interval (seconds)

    Returns
    -------

    tr,vr : array-like
        Regularized series
    """
    tr = np.arange(tnum[0], tnum[-1] + interval, interval)
    f = scipy.interpolate.interp1d(tnum, v, axis=0, assume_sorted=True, fill_value='extrapolate')
    vr = f(tr)
    return tr, vr


def nan_gap(t, h, interval=None):
    """ Ensure regular sampling for a time series and fill gaps with nans.
    NOTE: Not recommended if the series is intended for subsequent tidal analysis, e.g. the method doesn't work
    if base sampling changes, e.g. from sampling on 7 min of the hour to 27 min of the hour.

    interval : datetime.timedelta, optional
        Sampling interval. Determined from the time vector by default.
    """
    if interval is None:
        interval = find_si_robust(t)
    # work with float times: setdiff1d and intersect1d are much faster with floats than with datetimes
    tnum = numtime(t)
    interval_num = interval.total_seconds()
    tr, hr = regularize_num(tnum, h, interval_num)
    # fill gaps with nans
    # The below setdiff1d will not pick values at slightly irregular interval.
    # Solution is to fix those intervals or setdiff1d with tolerance.
    # tnan = np.setdiff1d(tr, tnum, assume_unique=True)
    # _, inv, _ = np.intersect1d(tr, tnan, return_indices=True)
    inv = setdiff_close(tr, tnum, tol=0.03*interval_num)  # 3% of the interval tolerance
    hr[inv] = np.nan
    return numtime_to_datetime(tr), hr


def setdiff_close(a, b, tol=1):
    """ Indices of elements in `a` that are not in `b` with tolerance """
    # Based on:
    # https://stackoverflow.com/questions/49303679/intersection-between-two-multi-dimensional-arrays-with-tolerance-numpy-pytho

    # Get closest distances for each pt in a
    # k=1 selects closest one neighbor
    dist = cKDTree(pkg_interp.atleast_2d0(b)).query(pkg_interp.atleast_2d0(a), k=1)[0]

    # Check the distances against the given tolerance value
    return dist > tol


def nan_gap2(t, h, interval=None, min_gap=None):
    """ Ensure regular sampling for a time series and fill gaps with nans.

    Interpolates to standard regular intervals, e.g. hourly on 00 of the hour, 15 min on 15 min of the hour etc.
    Recommended over nan_gap if the series is intended for subsequent tidal analysis, e.g. it works
    if base sampling changes, e.g. from sampling on 7 min of the hour to 27 min of the hour.

    interval : datetime.timedelta, optional
        Sampling interval. Determined from the time vector by default.

    Parameters
    ----------
    t : array_like of datetimes
    h : array_like
        Time series
    interval : datetime.timedelta; optional
        Sampling interval. Determined from the time vector by default.
    min_gap : datetime.timedelta; optional
        Shortest gap to fill with nans. Gaps shorter that that will be linearly interpolated over.
        Default is 1.1 * `interval`.

    Returns
    -------
    tr : array_like of datetimes
    hr : array_like
    """
    if interval is None:
        interval = find_si_robust(t)
    if min_gap is None:
        min_gap = interval * 1.1  # default in insert_nans_in_gaps
    tn, hn = insert_nans_in_gaps(t, h, interval=interval, tol=min_gap / interval)
    # round to the minute
    round_seconds = 60
    rounded_interval = datetime.timedelta(seconds=np.around(interval.total_seconds() / round_seconds) * round_seconds)
    tr = standard_time_vector(t, rounded_interval)
    f = scipy.interpolate.interp1d(numtime(tn), hn, axis=0, assume_sorted=True, fill_value='extrapolate')
    hr = f(numtime(tr))
    return tr, hr


def standard_time_vector(t, interval):
    """ Generate a time vector covering times in `t` with a specified interval.
    Rounds up to the nearest `interval` counting from datetime.datetime.min, i.e. datetime(1, 1, 1, 0, 0):

    Reference:
    https://stackoverflow.com/questions/32723150/rounding-up-to-nearest-30-minutes-in-python

    Parameters
    ----------
    t : array_like of datetimes
    interval : datetime.timedelta

    Returns
    -------
    array_like of datetimes at standard intervals
    """
    def ceil_dt(dt, delta):
        rounded_dt = dt + (datetime.datetime.min - dt) % delta
        # for timezone-aware datetime objects see the note at the link above
        localize = getattr(rounded_dt.tzinfo, 'localize', None)
        if localize:
            rounded_dt = localize(rounded_dt.replace(tzinfo=None), is_dst=bool(rounded_dt.dst()))
        return rounded_dt

    t0 = ceil_dt(t[0], interval)
    t1 = ceil_dt(t[-1], interval)
    return np.arange(t0, t1, interval).astype(datetime.datetime)


def insert_nans_in_gaps(t, h1, h2=None, interval=None, tol=1.1):
    """ Insert nan in each gap.
    Useful for further interpolation.

    Parameters
    ----------
    t : array_like; dtype=datetime
        Time vector.
    h1 : array_like
        Values.
    h2 : array_like (optional)
        Values. Use h1 & h2 for T&S or U&V, else use just h1 (SSH, T, ...)
    interval : datetime.timedelta; optional
        Sampling interval. Determined from the time vector by default.
    tol : float
        Tolerance for determining the gaps. Gap is a time interval which is greater than `interval` * `tol`.
    """
    if interval is None:
        interval = find_si_robust(t)
    if tol <= 1:
        raise RuntimeError('tolerance must be > 1')
    igap = np.diff(t) > interval * tol
    if np.any(igap):
        ta = np.append(t, t[:-1][igap] + interval)
        ha1 = np.append(h1, np.full(igap.sum(), np.nan))
        if h2 is not None:
            ha2 = np.append(h2, np.full(igap.sum(), np.nan))
        isort = np.argsort(ta)
        ts, hs1 = ta[isort], ha1[isort]
        if h2 is not None:
            hs2 = ha2[isort]
    else:
        ts, hs1, hs2 = t, h1, h2
    if h2 is None:
        return ts, hs1
    else:
        return ts, hs1, hs2


def find_si(t):
    """ Find sampling interval as the most frequent value in diff(t)
    
    Parameters
    ----------
    
    t : array_like
        Vector of numbers or datetime objects

    Returns
    -------
    
    Sampling interval, a number or datetime.timedelta
    """
    if not issorted(t):
        raise ValueError('Input must be sorted in ascending order.')
    return scipy.stats.mode(np.diff(t))[0][0]


def find_si_robust(t, return_unique_intervals=False):
    """ Robust determination of sampling interval
    
    Finds most frequent sampling interval after clustering of close values.
    
    Parameters
    ----------
    t : array_like
        Vector of datetime objects
    return_unique_intervals : bool, optional
        Return unique intervals between `t` values (in seconds). Close intervals are clustered.

    Returns
    -------
    Sampling interval, datetime.timedelta
    Unique intervals, seconds
    """
    
    if not issorted(t):
        raise ValueError('Input must be sorted in ascending order.')
    c = Counter(np.diff(t)) # count occurencies of each interval
    if len(c)==1:
        si_sec = next(iter(c.keys())).total_seconds() # Python 3
        si_cl = np.array([si_sec])
    else:
        counts = [k for k in c.values()]
        si = [k.total_seconds() for k in c.keys()] # sampling interval, decimal seconds
        # sort for subsequent clustering
        isort = np.argsort(si)
        si = np.array(si)[isort]
        counts = np.array(counts)[isort]
        # cluster very close si values (within 1% of most frequent si)
        # cl = scipy.cluster.hierarchy.fclusterdata(si[:,None],
        #                                           si[np.argmax(counts)]*0.01,
        #                                           criterion='distance')
        # cluster very close si values (within 8% of most frequent si); this deals with 1950 base year problem in nc files
        cl = scipy.cluster.hierarchy.fclusterdata(si[:, None],
                                                  si[np.argmax(counts)] * 0.08,
                                                  criterion='distance')
        counts_cl = accum_np(cl,counts) # sum counts within each cluster
        si_cl = accum_np(cl,si*counts)/counts_cl # weighted average of si within each cluster
        si_sec = si_cl[np.argmax(counts_cl)]
        
    # convert to timedelta
    # microseconds is the highest precision in timedelta, so round to 6 decimals
    si_sec = round(si_sec,6)
    sec = datetime.timedelta(seconds = int(si_sec // 1),
                             microseconds = int(round((si_sec % 1)*1e6)))
    if return_unique_intervals:
        return sec, si_cl
    else:
        return sec


def issorted(a):
    """ Check if vector is sorted in ascending order
    """
    return np.all(a[:-1] <= a[1:])


def accum_np(accmap, a, func=np.sum):
    """ Accumulate values in `a` according to the map in `accmap`

    Careful: This quick hack only works with contiguous accmaps,
    like 222111333, but not 1212323. Every change from one number to another
    will be seen as a new value. This avoids the slow sorting.

    Code from:
    https://mldesign.net/blog/2013/02/18/speedy-numpy-replacement-for-matlab-accumarray/
    """
    indices = np.where(np.ediff1d(accmap, to_begin=[1], to_end=[1]))[0]
    vals = np.zeros(len(indices) - 1)
    for i in range(len(indices) - 1):
        vals[i] = func(a[indices[i]:indices[i + 1]])
    return vals


def fill_small_gaps(t, v, max_gap=None):
    """ Fill small gaps (series of nans) in a time series by linear interpolation
    
    Parameters
    ----------
    
    t : (N,) array_like
        Vector of datetime objects
    v : (N,...) array_like
        Corresponding array of values with time along 1st dimension, 2D at most.
    max_gap : int, optional
        Longest gap to fill (number of data values). By default fills gaps of any length.

    Returns
    -------
    
    vi : array_like
        Series with filled small gaps, same shape as `v`
    """
    if v.ndim > 2:
        raise RuntimeError('Can work on 1D or 2D arrays only.')

    if v.ndim > 1:
        vi = v.copy()
        for k, v1 in enumerate(v.T):  # iterate over columns
            vi[:, k] = fill_small_gaps(t, v1, max_gap=max_gap)
    else:
        # fill all gaps
        ival = ~np.isnan(v)
        vi = np.interp(numtime(t), numtime(t[ival]), v[ival])
        #TODO What about leading or trailing NaNs?

        if max_gap is not None:
            # fill large gaps back with nans
            ist,ien = continuous_chunks(np.isnan(v))
            for k in range(len(ist)):
                if ien[k]-ist[k] > max_gap:
                    vi[ist[k]:ien[k]] = np.nan # fill with nans
    
    return vi


def continuous_chunks(v):
    """ Find continuous intervals of 1's or True's in a vector of 1/0 or True/False
    
    Parameters
    ----------
    
    v : array_like
        Vector of 0/1 or True/False.
    
    Returns
    -------
    
    ist,ien : array_like, int
        Start and end indices of continuous intervals.
        
    Examples
    --------
    Find indices of gaps in a regularly spaced series with gaps filled with nans
    
    ist,ien = continuous_chunks(np.isnan(x))
    
    Find indices of continuous intervals of good data
    
    ist,ien = continuous_chunks(~np.isnan(x))
    
    """
    v = v.astype(int) # ensure 0/1 values 
    dn = np.diff(v)
    ist = np.where(dn==1)[0] + 1 # chunk start (zero-based)
    ien = np.where(dn==-1)[0] + 1 # chunk end (zero-based, exclusive)
    
    if v[0]==1:
        ist = np.r_[0,ist]
    if v[-1]==1:
        ien = np.r_[ien,len(v)]
    
    return ist,ien


def filterdata(x, fco, fs=1, rolloffband=None, attenuation=1e-4, showinfo=False):
    """Low-pass filtering

    Low-pass Kaiser-windowed FIR filter is designed with optimal parameters
    estimated with kaiserord. Similar to Matlab FIR filtering approach.

    Parameters
    ----------

    x : array_like
        Series to filter. Works on columns of `x` in case it is a 2D array.
    fco : float
        Cutoff frequency (in same units as `fs`)
    fs : float, optional
        Sampling frequency, default is 1.
    rolloffband : float, optional
        FIR filter parameter: approximate 1/2 width of rolloff in time
        units (inverse of `fs` units)
    attenuation : float, optional
        FIR filter parameter: maximum ripple in pass band and
        minimum attenuation in stop band. Attenuation is relative energy
        suppression level (not in decibels). Default is 1e-4, i.e. the energy
        in stop band will be suppressed at least 1e-4 times of the original level.
    showinfo : bool, optional
        True to display window legth and beta, False for silent run.
        Default is False.

    Returns
    -------

    y : array_like
        Filtered series.

    Examples
    --------
    Filter a series with 30-min sampling, cutoff 1/20 cpd,
    rolloff from 1/18 to 1/22 cpd, energy attenuation in stop band 10^-4::

      y = filterdata(x, 1/20, fs=48, rolloffband=2, attenuation=1e-4)

    """
    if rolloffband is None:
        rolloffband = 1 / fco / 12  # half-rolloff, e.g. 2/24 for daily resampling with 30-hr cutoff

    width = 1 / (1 / fco - rolloffband) - 1 / (1 / fco + rolloffband)
    #    width = 2 * (fco - 1/(1/fco+rolloffband)) # slightly narrower band
    winlen, beta = scipy.signal.kaiserord(-20 * np.log10(attenuation), width / (fs / 2))
    winlen += -(winlen % 2) + 1  # ensure odd length
    if showinfo:
        print('Filtering with window length: {}, beta: {}'.format(winlen, beta))
    b = scipy.signal.firwin(winlen, fco, window=('kaiser', beta), fs=fs)

    if x.ndim == 1:
        x = x[:, None]
        x_is_vector = True
    else:
        x_is_vector = False

    ncol = x.shape[1]

    y = np.zeros(x.shape)
    for k in range(ncol):  # convolution for each column
        y[:, k] = convolve_nan(x[:, k], b)

    if x_is_vector:
        y = y[:, 0]  # preserve shape of the input

    return y


def convolve_nan(x, b):
    """ Convolve a vector (with NaNs) with a window

    Parameters
    ----------
    x : array_like
        Vector of values
    b : array_like
        Window coefficients, must be odd length

    Returns
    -------
    y : array_like
        Convolved vector of same size as x
    """
    nx, nb = len(x), len(b)
    # pad with nans on both sides
    npad = int((nb - 1) / 2)
    pad = np.nan * np.ones(npad)
    xp = np.r_[pad, x, pad]
    # keep track of nans and fill zeros for np.sum to work
    inv = np.isnan(xp)
    ival = ~inv
    xp[~ival] = 0
    y = np.zeros_like(x)
    # do convolution with window scaling
    for k in range(nx):
        xk = xp[k:k + nb]
        ivalk = ival[k:k + nb]
        y[k] = np.sum(xk * b) / np.sum(b[ivalk])
    # restore original nans
    y[inv[npad:npad + nx]] = np.nan
    return y


def resample_filt(t, v, target_si_days, midflag=False, fill_gaps=True):
    """ Generic resampling routine. Downsamples after low-pass filtering.
    Works on columns of b.

    Parameters
    ----------
    t : (N,) array_like
        Time vector.
    v : (N,...) array_like
        Corresponding array of values with time along 1st dimension, 2D at most.
    target_si_days : float
        Required sampling interval (days).
    midflag : bool
        Shift sampling to the middle of intreval; can be useful when resampling to daily or longer intervals.

    Returns
    -------

    """
    si = find_si_robust(t)
    if fill_gaps:
        tr, vr = regularize(t, v, si)  # fill sampling gaps linearly, can still contain nans
        vf = fill_small_gaps(tr, vr)  # fill all nans linearly
    else:
        tr = t
        vf = v.copy()

    # % % remove leading or trailing NaNs
    # % %%% brute force: remove all values with NaNs in ANY column
    # % anyinv = any(isnan(b'));
    # % t(anyinv) = [];
    # % b(anyinv,:) = [];
    # % inv(anyinv,:) = [];
    # %
    # % si = find_si(t);
    # % si = round(si*24*60*60)/24/60/60; % round to the second
    # %
    # % bf(inv) = NaN;

    fs = 1 / (si.total_seconds() / 24 / 60 / 60)  # sampling freq (cpd)
    fco = 1 / target_si_days / 1.25  # fco corresponding to target sampling interval +25%
    vff = filterdata(vf, fco, fs=fs)

    target_delta = datetime.timedelta(days=target_si_days)
    td = np.arange(tr[0], tr[-1] + target_delta, target_delta).astype(datetime.datetime)

    if midflag and target_si_days >= 1:
        # for resampling to daily or longer time steps, move sampling to the middle
        # of the interval, e.g. sample daily data at noon
        td += target_delta / 2
        td = np.r_[td[0] - target_delta, td]  # ensure 1st value is in
        td = td[(td >= t[0]) & (td <= t[-1])]  # all time stamps covered by the original timeline

    # remove times above the interpolation range
    td = td[td <= tr[-1]]
    f = scipy.interpolate.interp1d(numtime(tr), vff, axis=0, assume_sorted=True)
    vd = f(numtime(td))

    return td, vd


def spikes(v, param=10, iterate=True, verbose=True):
    """ Detect spikes in a vector v by a threshold relative to standard deviation

    It is advised to pass a high-pass filetered series to this functon.

    Parameters
    ----------

    v : array_like
        Vector to despike (regularly sampled)
    param : float or int, optional
        Spike threshold parameter, multiples of standard deviation (default is 10)
    iterate : bool, optional
        Use iterative approach: run spike detection until no spikes detected in
        the de-spiked series (default is True)
    verbose : bool, optional
        Print info on spikes detection. Default is True

    Returns
    -------

    inv : array_like
        Logical indices of spikes
    """

    sdev = np.nanstd(v, axis=0)
    smean = np.nanmean(v, axis=0)
    with np.errstate(invalid='ignore'):
        inv = np.abs(v - smean) > sdev * param  # spike indices at this iteration
    if verbose:
        print('spikes: First iteration')
        print('   Max deviation = ', np.nanmax(np.abs(v - smean)) / sdev, ' of standard deviation')
        print('   Spikes detected: ', np.sum(inv))
    if iterate:
        vn = v.copy()
        vn[inv] = np.nan
        invk = inv.copy()
        while np.any(invk):
            sdev = np.nanstd(vn, axis=0)
            smean = np.nanmean(vn, axis=0)
            with np.errstate(invalid='ignore'):
                invk = np.abs(vn - smean) > sdev * param  # spike indices at this iteration
            vn[invk] = np.nan
            inv = np.logical_or(inv, invk)  # add newly detected spikes
            if verbose:
                print('spikes: Iterating...')
                print('   Max deviation = ', np.nanmax(np.abs(vn - smean)) / sdev, ' of standard deviation')
                print('   Spikes detected: ', np.sum(invk))
    return inv


def zero_crossings(x, y):
    """ Zero-crossings of a linear segmented curve y = f(x).

    Parameters
    ----------

    x,y : array-like
        1D arrays, x-coordinate (can be array of datetime) and corresponding y values

    Returns
    -------

    x0 : array-like
        Coordinates of zero-crossings
    sign : array-like
        Corresponding direction of change:
        +1 for change of sign from negative to positive,
        -1 for change of sign from positive to negative,
        0 if curve touches, but does not cross zero
    """
    difsign = np.diff(np.sign(y))
    difsign[np.isnan(difsign)] = 0  # exclude changes to and from nan
    i1 = np.where(difsign != 0)[0]  # points before zero-crossings
    i2 = i1 + 1  # points after zero-crossings
    x0 = x[i1] + (x[i2] - x[i1]) * (-y[i1]) / (y[i2] - y[i1])  # linear interpolation
    sign = np.sign(y[i2])
    sign[np.sign(y[i1]) == 0] = 0  # needed if 1st point is zero or if curve touches zero
    x0, iu = np.unique(x0, return_index=True)  # needed if curve touches zero
    sign = sign[iu]  # also keep sign only for unique points
    return x0, sign

