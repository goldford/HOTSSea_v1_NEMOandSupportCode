# created 2023-08-17 by GO
# keep helper code functions etc out of notebooks
from datetime import datetime, timedelta
import numpy as np
import calendar
from stats_GO import willmott1981
import skill_metrics as sm


def seconds_in_one_month(year, month):
    _, num_days = calendar.monthrange(year, month)
    return num_days * 24 * 60 * 60


# generic mvg avg filter using user-defined months and / or days
def apply_mvgavg_filter(time_sorted, var_sorted, mos=1, dys=30, use_nanmean=True):
    three_months_in_seconds = mos * dys * 24 * 60 * 60

    result_time = []
    result_var = []

    start_time = datetime.fromtimestamp(time_sorted[0])
    end_time = datetime.fromtimestamp(time_sorted[-1])

    current_time = start_time
    while current_time <= end_time:
        # Calculate the end of the 4-month period
        end_of_period = current_time + timedelta(seconds=three_months_in_seconds)

        # Create a time mask for the current 4-month period
        time_mask = (time_sorted >= current_time.timestamp()) & (time_sorted < end_of_period.timestamp())

        # Calculate the average var concentration for the current 4-month period
        if use_nanmean:
            average_var = np.nanmean(var_sorted[time_mask])
        else:
            average_var = np.mean(var_sorted[time_mask])

        result_time.append(current_time.timestamp())
        result_var.append(average_var)

        current_time += timedelta(seconds=three_months_in_seconds)

    # Convert the result_time back to datetime objects
    result_time = [datetime.fromtimestamp(ts) for ts in result_time]

    # Convert the result arrays to numpy arrays
    result_time = np.array(result_time)
    result_var = np.array(result_var)
    return result_time, result_var


# get stats for taylor fig from SST scores from buoys
# bu - buoy
# mv_avg - in days
# SSTbuoy_scores - from pyap (loaded pickle file output from analyse.py)
# run_sname - alt name for model run
# use_nanmean - not used!
def get_taylor_stats_SST(stn, mv_avg, SST_scores, run_sname, use_nanmean):
    # if 1 day then just use the precalculated pyap stats

    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        ccoef = scores['pearson']
        wss = scores['skill1981']

        # for target
        rmsd = scores['rmse']
        bias = scores['bias']
        nrmsd = rmsd / stdev_obs
    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o['temperature'].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m['temperature'].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        stats = sm.taylor_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, 'data')
        sdev = stats['sdev'][1]
        crmsd = stats['crmsd'][1]
        ccoef = stats['ccoef'][1]
        wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        mod_norm_stdev = stdev_mod / stdev_obs

        ncrmsd = crmsd / stdev_obs

    return [stn, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss]


def get_taylor_stats_LH(stn, mv_avg, SST_scores, run_sname, use_nanmean, var_ts):
    # get stats for taylor fig from scores from LH
    # created by HGo 2023-08-18 - only diff from above is
    #    there are both salin and temp data
    # stn - lighthouse
    # mv_avg - in days
    # SSTbuoy_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # use_nanmean - not used!

    # if 1 day then just use the precalculated pyap stats

    if var_ts == "T":
        var_ = 'temperature'
    elif var_ts == "S":
        var_ = 'salinity'
    else:
        print("issue with var passed to get_taylor_stats")

    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores'][var_]
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        ccoef = scores['pearson']
        wss = scores['skill1981']

    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o[var_].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m[var_].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        stats = sm.taylor_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, 'data')
        sdev = stats['sdev'][1]
        crmsd = stats['crmsd'][1]
        ccoef = stats['ccoef'][1]
        wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        mod_norm_stdev = stdev_mod / stdev_obs

        ncrmsd = crmsd / stdev_obs

    return [stn, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss]


def get_target_stats_SST(stn, mv_avg, SST_scores, run_sname, use_nanmean=True, augment_rmsd=True):
    # get stats for target fig from SST scores
    # stn - buoy
    # mv_avg - in days
    # SST_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # augment_rmse - multiply RMSD by sign of mod_stdev - obs
    # use_nanmean - not used!
    # if 1 day then just use the precalculated pyap stats
    if mv_avg == 1:
        scores = SST_scores[run_sname][stn]['filt_interp_data']['scores']
        bias = scores['bias']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        # mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        rmsd = scores['rmse']
        nrmsd = rmsd / stdev_obs



    else:
        filt_interp_data_o = SST_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = SST_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o['temperature'].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m['temperature'].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        # returns (optionally normalized)
        # {'bias': bias, 'crmsd': crmsd, 'rmsd': rmsd}
        stats = sm.target_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, norm=True)
        bias = stats['bias']
        ncrmsd = stats['crmsd']
        nrmsd = stats['rmsd']

        # sdev = stats['sdev'][1]
        # crmsd = stats['crmsd'][1]
        # # ccoef = stats['ccoef'][1]
        # # wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        # # mod_norm_stdev = stdev_mod / stdev_obs
        # ncrmsd = crmsd / stdev_obs

    if augment_rmsd:
        if (stdev_mod - stdev_obs) < 0:
            ncrmsd = ncrmsd * -1

    return [stn, bias, ncrmsd, nrmsd]


def get_target_stats_LH(stn, mv_avg, LH_scores, var_ts, run_sname, use_nanmean=True, augment_rmsd=True):
    # get stats for target fig from LH scores
    # stn - buoy
    # mv_avg - in days
    # LH_scores - from pyap (loaded pickle file output from analyse.py)
    # run_sname - alt name for model run
    # augment_rmse - multiply RMSD by sign of mod_stdev - obs
    # use_nanmean - not used!

    if var_ts == "T":
        var_ = 'temperature'
    elif var_ts == "S":
        var_ = 'salinity'
    else:
        print("issue with var passed to get_taylor_stats")

    if mv_avg == 1:
        scores = LH_scores[run_sname][stn]['filt_interp_data']['scores'][var_]
        bias = scores['bias']
        stdev_obs = scores['stdev_obs']
        stdev_mod = scores['stdev_mod']
        # mod_norm_stdev = stdev_mod / stdev_obs
        crmsd = scores['crmse']
        ncrmsd = crmsd / stdev_obs
        rmsd = scores['rmse']
        nrmsd = rmsd / stdev_obs

    else:
        filt_interp_data_o = LH_scores['obs'][stn]['filt_interp_data']
        filt_interp_data_m = LH_scores[run_sname][stn]['filt_interp_data']
        # get original time series, calculate moving average
        time_obs = np.asarray(filt_interp_data_o['time'])
        t_obs = np.asarray(filt_interp_data_o[var_].astype('f'))
        time_mod = np.asarray(filt_interp_data_m['time'])
        t_mod = np.asarray(filt_interp_data_m[var_].astype('f'))

        # format to datetime
        time_obs_numeric = np.array([dt.timestamp() for dt in time_obs])
        time_mod_numeric = np.array([dt.timestamp() for dt in time_mod])

        # sort
        sort_indices = np.argsort(time_obs_numeric)
        time_sort_o = time_obs_numeric[sort_indices]
        t_sort_o = t_obs[sort_indices]

        sort_indices = np.argsort(time_mod_numeric)
        time_sort_m = time_mod_numeric[sort_indices]
        t_sort_m = t_mod[sort_indices]

        time_obs_mvgavg, t_obs_mvgavg = apply_mvgavg_filter(time_sort_o, t_sort_o, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)
        time_mod_mvgavg, t_mod_mvgavg = apply_mvgavg_filter(time_sort_m, t_sort_m, mos=1,
                                                            dys=mv_avg, use_nanmean=use_nanmean)

        # stats = {'ccoef': ccoef, 'crmsd': crmsd, 'sdev': sdev}
        # remove nans to avoid crash
        t_mod_mvgavg_nonan = t_mod_mvgavg[~np.isnan(t_mod_mvgavg) & ~np.isnan(t_obs_mvgavg)]
        t_obs_mvgavg_nonan = t_obs_mvgavg[~np.isnan(t_obs_mvgavg) & ~np.isnan(t_mod_mvgavg)]

        # returns (optionally normalized)
        # {'bias': bias, 'crmsd': crmsd, 'rmsd': rmsd}
        stats = sm.target_statistics(t_mod_mvgavg_nonan, t_obs_mvgavg_nonan, norm=True)
        bias = stats['bias']
        ncrmsd = stats['crmsd']
        nrmsd = stats['rmsd']

        # sdev = stats['sdev'][1]
        # crmsd = stats['crmsd'][1]
        # # ccoef = stats['ccoef'][1]
        # # wss = willmott1981(t_obs_mvgavg, t_mod_mvgavg)
        stdev_obs = np.nanstd(t_obs_mvgavg)
        stdev_mod = np.nanstd(t_mod_mvgavg)
        # # mod_norm_stdev = stdev_mod / stdev_obs
        # ncrmsd = crmsd / stdev_obs

    if augment_rmsd:
        if (stdev_mod - stdev_obs) < 0:
            ncrmsd = ncrmsd * -1

    return [stn, bias, ncrmsd, nrmsd]
