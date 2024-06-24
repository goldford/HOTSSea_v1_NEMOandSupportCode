import numpy as np
import scipy as sp
import scipy.stats

# changelog: 
# GO 20230919 - added mean val calc over deps that accounts for stretched z levs

def residual_stats(obs, mod):
    """
    Error statistics for scalar data (water level, temperature, salinity, ...)
    """

    if isinstance(obs[0], complex):
        # temporary redirect here until we update all velocity analysis code to call residual_stats_vector directly
        return residual_stats_vector(obs, mod)

    if np.all(np.isnan(obs)) or np.all(np.isnan(mod)):
        return {'skill1981': np.nan, 'mean_obs': np.nan, 'mean_mod':np.nan, 'bias': np.nan, 'crmse': np.nan, 'rmse': np.nan,
                'gamma2': np.nan, 'pearson': np.nan, 'mae': np.nan, 'mad': np.nan,'stdev_obs': np.nan, 'stdev_mod':np.nan}

    scores = {}
    scores['skill1981'] = Willmott1981(obs, mod)
    scores['mean_obs'] = np.nanmean(obs)
    scores['mean_mod'] = np.nanmean(mod)
    scores['bias'] = np.nanmean(mod - obs)  # why not use the bias() function below?
    scores['crmse'] = np.nanstd(mod - obs)
    scores['gamma2'] = np.nanvar(mod - obs) / np.nanvar(obs)
    scores['rmse'] = np.nanmean((mod - obs) ** 2) ** 0.5
    # Mean Absolute Error:
    scores['mae'] = np.nanmean(np.abs(mod-obs))
    # Mean Absolute Deviation:
    scores['mad'] = np.nanmean(np.abs(mod-obs - np.nanmean(mod-obs)))

    i = ~np.isnan(mod-obs)
    if sum(i) < 2: # pearsonr is ill-defined for less than 2 values
        scores['pearson'] = np.nan
    else:
        scores['pearson'], _ = sp.stats.pearsonr(mod[i], obs[i])

    # GO added May 2023 for Taylor diag
    scores['stdev_obs'] = np.nanstd(obs)
    scores['stdev_mod'] = np.nanstd(mod)

    return scores


def residual_stats_vector(obs, mod):
    """
    Error statistics for complex data (velocities) where velocity = u + 1j*v
    """

    if np.all(np.isnan(obs)) or np.all(np.isnan(mod)):
        return {'mean_obs': np.nan, 'mean_mod': np.nan, 'bias': np.nan, 'crmse': np.nan, 'rmse': np.nan, 'gamma2': np.nan,
                'vector_correlation': np.nan, 'vector_correlation_magnitude': np.nan,
                'vector_correlation_radians': np.nan, 'vector_correlation_degrees': np.nan}
    scores = {}

    scores['mean_obs'] = np.nanmean(obs)
    scores['mean_mod'] = np.nanmean(mod)
    scores['bias'] = np.nanmean(mod - obs)  # why not use the bias() function below?
    scores['crmse'] = np.nanstd(mod - obs)
    # Cummins and Thupaki 2018 Eqn 6 - vector rms error
    scores['rmse'] = np.nanmean( np.abs(mod - obs) ** 2) ** 0.5
    scores['gamma2'] = np.nanvar(mod - obs) / np.nanvar(obs)

    # vector correlation
    vcc = Rohrs2015_vector_correlation(obs, mod)
    scores['vector_correlation'] = vcc
    scores['vector_correlation_magnitude'] = np.absolute(vcc)
    scores['vector_correlation_radians'] = np.angle(vcc, deg=False)
    scores['vector_correlation_degrees'] = np.angle(vcc, deg=True)

    return scores


def circular_stats (obs,mod):
    """Used for calculating statistics on things like direction, which vary from 0 to 2pi. """
    # CAUTUION: THIS FUNCTION IS NOT FUNCTIONAL YET
    if np.all(np.isnan(obs)) or np.all(np.isnan(mod)):
        scores = {'skill1981': np.nan, 'bias': np.nan, 'crmse': np.nan, 'rmse': np.nan,
                  'gamma2': np.nan, 'vector_r': np.nan, 'pearson': np.nan, 'mae': np.nan, 'mad': np.nan}
    else:
        i = ~np.isnan(mod-obs)
        scores = {}
        scores['bias'] = sp.stats.circmean(i)
        scores['crmse'] = sp.stats.circstd(i)
        scores['var']  = sp.stats.circvar(i)
        scores['mae']  = sp.stats.circmean(np.abs(i))
        scores['mad']  = sp.stats.circmean(np.abs(i - np.nanmean(i) ) )  # should this subtract sp.stats.circmean(i) instead of np.nanmean(i)?
        scores['rmse'] = (sp.stats.circmean(i**2)) **0.5
        scores['pearson'] = np.nan
        scores['vector_r'] = np.nan

        # Correlation from page 8 of https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Correlation.pdf
        omean = sp.stats.circmean(obs)
        mmean = sp.stats.circmean(mod)
        numer = np.sum(np.sin(obs-omean)*np.sin(mod-mmean) )
        denom = np.sum(np.sin(obs-omean)**2 ) * np.sum(np.sin(mod-mmean)**2)
        scores['correlation'] = numer / np.sqrt(denom)

    return scores




# Following two metrics are from Cummins and Thupaki 2018 doi:10.1016/j.csr.2017.10.007
# Ellipse errors are their own, tidal error is used for water level and cited as the starting point
# Total, counterotating, and relative errors are all included here
def ellipse_error(obs, mod, constituents):
    """Calculates the error of the modelled ellipse relative to the observed.  
    All four tidal parameters are included in the calculation, making it a comprehensive 
    evaluation of the skill of the model.  
    Errors from each of the two counterrotating ellipses are also calculated, as is the relative
    error (ie the absolute error scaled by a measurement of the observed ellipse).
    Citation for all quantities: Cummins and Thupaki 2018, doi:10.1016/j.csr.2017.10.007"""

    # obs and mod input parameters are tuples of ellipse parameters for each of the constituents in "constituents"
    # obs = (omaj, omin, oinc, opha) and similarly for mod
    # Notation from Cummins and Thupaki 2018 has semi major/minor axes as A/B, inclination as theta, phase lage as g.
    # Unpack the tuple using their notation for clarity of comparison with published work.

    Du, Durel, Dccw, Dcw = {}, {}, {}, {}

    for con in constituents:
        # Unpack the tuples
        oA, oB, otheta, og = obs[con]
        mA, mB, mtheta, mg = mod[con]

        otheta = np.radians(otheta)
        og = np.radians(og)
        mtheta = np.radians(mtheta)
        mg = np.radians(mg)

        # Cummins and Thupaki 2018 Eqn 8
        Du[con] = np.sqrt(
            0.5*(oA**2 + oB**2 + mA**2 + mB**2)
            - np.cos(og-mg) * np.cos(otheta-mtheta) * (oA*mA+oB*mB)
            - np.sin(og-mg) * np.sin(otheta-mtheta) * (oA*mB+mA*oB)
        )

        # Cummins and Thupaki 2018 Eqn 13
        denom = np.sqrt((oA**2 + oB**2)/2.0)
        Durel[con] = Du[con] / denom

        # Cummins and Thupaki 2018 Eqn 11 (as well as 10)
        oap = 0.5*(oA+oB);    oam = 0.5*(oA-oB);    oep = otheta - og;    oem = otheta+ og
        map = 0.5*(mA+mB);    mam = 0.5*(mA-mB);    mep = mtheta - mg;    mem = mtheta+ mg
        Dccw [con] = np.sqrt ( oap**2+map**2 - 2.0*oap*map*np.cos(oep-mep)  )    
        Dcw  [con] = np.sqrt ( oam**2+mam**2 - 2.0*oam*mam*np.cos(oem-mem)  )    

    return Du, Durel, Dccw, Dcw


# Reference in Cummins and Thupaki 2018
def tidal_error(obs, mod, constituents):
    """ Calculates the total tidal error for water level, including effect of both amplitude and phase. """
    delta = {}
    for con in constituents:
        # Unpack the tuples, using Cummins and Thupaki 2018 notation
        oh, ophi, _, _ = obs[con]
        mh, mphi, _, _ = mod[con]
        # Cummins and Thupaki 2018 Eqn 1, ref'd from Cummins and Oey 1997
        delta[con] = np.sqrt(0.5*(oh**2 + mh**2) - oh * mh * np.cos(np.radians(ophi) - np.radians(mphi)))
    return delta


def tidal_error_common(obs_tide_result, mod_tide_result):
    """ Tidal error for constituents common to two tidal analysis results
    according to Cummins and Thupaki 2018 Eqn 1, ref'd from Cummins and Oey 1997.
    Parameters
    ----------
    obs_tide_result : dict
        Output from t_tide for the 1st series, e.g. observations
    mod_tide_result : dict
        Output from t_tide for the 2nd series, e.g. model
    Returns
    -------
    Dictionary with constituent names as keys and tidal errors as values
    """
    _, io, im = np.intersect1d(obs_tide_result['fu'], mod_tide_result['fu'], return_indices=True)
    const_names = obs_tide_result['nameu'][io]
    oh = obs_tide_result['tidecon'][io, 0]
    ophi = obs_tide_result['tidecon'][io, 2]
    mh = mod_tide_result['tidecon'][im, 0]
    mphi = mod_tide_result['tidecon'][im, 2]
    errs = np.sqrt(0.5 * (oh ** 2 + mh ** 2) - oh * mh * np.cos(np.radians(ophi) - np.radians(mphi)))
    return dict(zip(const_names, errs))


# Confidence interval calculation.  Informed by Sasha's surge evaluation package
def stat_with_conf_interval(stat, obs, mod, axis=0):
    """ Model error statistics with 95% confidence interval estimate """
    if stat not in globals():
        raise RuntimeError(f'Function {stat} not found in pkg_statistics')
    fun = globals()[stat]
    s = fun(obs, mod, axis=axis)
    rng = np.random.default_rng()  # can also be left at default in bootstrap(), not much difference
    # 95% confidence by default
    res = scipy.stats.bootstrap((obs, mod), fun, axis=axis, paired=True, random_state=rng)
    return s, res.confidence_interval.low, res.confidence_interval.high

# ECCC stats calcs are at
# https://gitlab.science.gc.ca/olh001/surge_validation/blob/4f7cee04adb023af82b3c06bb0ca93a21bea2267/src/surge_validation/verification_stats/calc_stats_with_obs.py#L51


# Individual statistics functions
# Currently all for real numbers only -- may need generalization for complex numbers.
# All functions should accept two arrays and the axis for vectorized calculations.
# This ensures compatibility with conf_interval().
def gamma2(obs, mod, axis=0):
    """Calculate the gamma^2 = var(mod-obs) / var(obs) for ab obs-mod pair.
    Depending on how the obs and mod are sliced, this can return a function of time or ensemble."""
    return np.nanvar(np.asarray(mod)-np.asarray(obs), axis=axis) / np.nanvar(np.asarray(obs), axis=axis)


def rmse(obs, mod, axis=0):
    """ Calculate RMSE = sqrt (1/n sum (mod-obs)^2 ) for an obs-mod pair. """
    return np.nanmean((np.asarray(mod) - np.asarray(obs)) ** 2, axis=axis) ** 0.5


def bias(obs, mod, axis=0):
    return np.nanmean(np.asarray(mod), axis=axis) - np.nanmean(np.asarray(obs), axis=axis)


def crmse(obs, mod, axis=0):
    return np.nanstd(np.asarray(mod) - np.asarray(obs), axis=axis)

def pearson (obs, mod):
    i = ~np.isnan(np.asarray(mod)-np.asarray(obs))
    return sp.stats.pearsonr(np.asarray(mod[i]), np.asarray(obs[i]))   # no nanpearson


def mae(obs, mod, axis=0):
    return np.nanmean(np.abs(np.asarray(mod) - np.asarray(obs)), axis=axis)
    
# GO added 20230612
def nanmean_simple(d):
    return np.nanmean(d)

def nanstd_simple(d):
    return np.nanstd(d)


# added by GO 20230920, cross checked
# todo: this is meant only for CTD data analysis (only for 't' depths)
def nanmean_vvl(d, mDep, e3t):
    # d - mod or obs data interpolated to model depths (assumed 't')
    # mDep - model deps (truncated to d shape)
    # e3t0 - depth weights or spans vertically in metres, truncated to shape of d
    
    range = mDep[-1] - mDep[0]
    sum_e3t = np.nansum(e3t)
    if sum_e3t == 0:
      return np.nan
    weights = e3t / sum_e3t
    weighted_mean = np.nansum((d * weights),axis=0)
     
    
    return weighted_mean

# added by GO 20230920, cross checked
# todo: this is meant only for CTD data analysis (only for 't' depths)
def nanstd_vvl(d, mDep, e3t):
    # Calculate the weighted mean first
    range = mDep[-1] - mDep[0]
    sum_e3t = np.nansum(e3t)
    if sum_e3t == 0:
        return np.nan
    weights = e3t / sum_e3t
    weighted_mean = np.nansum((d * weights), axis=0)
    if np.isnan(weighted_mean).any():
        return np.nan
    # Calculate the weighted variance
    weighted_variance = np.nansum(weights * (d - weighted_mean)**2, axis=0)
    if np.isnan(weighted_variance).any():
        return np.nan
    # Calculate the weighted standard deviation
    weighted_std = np.sqrt(weighted_variance)
    
    return weighted_std


def Willmott1981(obs, mod, axis=0):
    num = np.nansum((mod - obs) ** 2, axis=axis)
    obs_mean = np.nanmean(obs, axis=axis)
    dM = np.abs(mod - obs_mean)
    dO = np.abs(obs - obs_mean)
    den = np.nansum((dM + dO) ** 2, axis=axis)
    if den == 0:
        return np.nan
    else:
        return np.max([0, 1 - num / den])


def Willmott1985(obs, mod, axis=0):
    num = np.nansum(np.abs(mod - obs), axis=axis)
    obs_mean = np.nanmean(obs, axis=axis)
    dM = np.abs(mod - obs_mean)
    dO = np.abs(obs - obs_mean)
    den = np.nansum((dM + dO), axis=axis)
    if den == 0:
        return np.nan
    else:
        return np.max([0, 1 - num / den])


def Willmott2012(obs, mod, axis=0):
    c = 2
    num = np.nansum(np.abs(mod - obs), axis=axis)
    dO = np.nansum(np.abs(obs - np.nanmean(obs, axis=axis)), axis=axis)
    if num <= c*dO:
        return 1 - (num/c*dO)
    else:
        return (c*dO/num) - 1


def Kundu1976_vector_correlation(obs, mod):
    """
    Equation 3.4 from Kundu 1976
    Reference:
      Kundu, P. K. (1976). Ekman Veering Observed near the Ocean Bottom, Journal of Physical Oceanography, 6(2), 238-242.
      https://journals.ametsoc.org/view/journals/phoc/6/2/1520-0485_1976_006_0238_evonto_2_0_co_2.xml
    """
    i = ~np.isnan(mod-obs)
    if sum(i) < 2: # correlation is ill-defined for less than 2 values
        return np.nan + 1j*np.nan
    else:
        rnum = np.nanmean(np.conj(mod) * obs)
        rden = np.nanmean(np.conj(mod) * mod) * np.nanmean(np.conj(obs) * obs)
        vc = rnum / np.sqrt(rden)
        # Now we convert from counterclowise-positive to clockwise-positive
        vc1 = np.conj(vc)
        return vc1


def Rohrs2015_vector_correlation(obs, mod):
    """
    Equation 6 from Röhrs 2015 -- same as Kundu, 1976 except we remove the means.
    Reference:
     Röhrs, J., and Christensen, K. H. (2015), Drift in the uppermost part of the ocean,
     Geophys. Res. Lett., 42, 10,349– 10,356, doi:10.1002/2015GL066733.
    """
    return Kundu1976_vector_correlation(obs - np.nanmean(obs), mod - np.nanmean(mod))
