
# Created Jan 2024 by G Oldford
# purpose: generate trend statistics for model and obs data from the Nanoose Station time series
# input:
#    meshm_p     - string; path to mesh mask from NEMO which contains depths and coords
#    meshm_f     - string; mesh mask file name
#    dat_f       - string; netcdf data file of observations, assumed biweekly and spatially averaged (dims: depth, time)
#    depth_min   - float; minimum depth to use (during depth integration)
#    depth_max   - float; maximum depth to use (during depth integration)
#    var         - string; the variable being analysed (e.g., "temperature") - only tested with temp for now
#    resolution  - float; the precision of the measurement (sig digits) for ties calculated during some stats calcs
#    remove_nans - True/False;
#    time_inc    - string; sub-annual time increment for temporal avg ("biweekly", "monthly", "seasonal", "annual")
#    alpha_CI    - the alpha confidence interval threshold (e.g., 0.05 for 95%)
#    nan_mask_in - in case we wish to filter out values using a mask from a different dataset
# returns:
#    slopes
#    summary_stats_all
#    DTS_CI_out        - results of trend test using decorr. time scale
#                      - list of dict elements; {season: {slope, inter, r_val, p_val, std_err, slope_lcl, slope_ucl}}
#    mk_3pw_out        - results of '3PW' trend test using Mann-Kendall, see https://doi.org/10.5194/amt-13-6945-2020
#                      - list of dict elements; {P= max(P_PW, P_TFPW_Y);'ss' (float): statistical significance:
#                          alpha_MK if the test is ss at the alpha confidence level. Defaults to 95.
#                          0 if the test is not ss at the alpha_MK confidence level.
#                          -1 if the test is a TFPW_Y false positive at alpha_MK confidence level
#                          -2 if the test is a PW false positive at alpha_MK confidence level
#                          'slope' (float): Sen's slope in units/y
#                           'ucl' (float): upper confidence level in units/y
#                           'lcl' (float): lower confidence level in units/y
# notes:
#    - the modelled data inputs may or may not be 'blind' (nans where obs has gaps)
#    - for time averaging, the finest resolution assumed present in the data is biweekly
#    - to do: bootstrap, make below a function that can be looped, auto-select trend test after stats done,
#             introduce some kind of harmonic analysis or ESD
#    - the nan_mask system in place implicitly assumes that input files have nans and the same
#      shape as each other! if this is called iteratively...


import os
import numpy as np
import xarray as xr
from GO_tools import (get_meshmask, depth_int, do_fft, do_seasonal_avg,
                      do_ss_seasons, do_ss, do_lr_seasons, prep_data_seas, do_summary_stats,
                      dst_trend_sig, mk3pw_trend_sig)
import datetime as dt
from tabulate import tabulate


def do_prelim_trend_analysis(dat_davg, dat_pytime, dat_numerictime,
                             dat_seas, dat_pytime_seas, dat_timenumeric_seas,
                             seasons_per_year, dat_ts, resolution, alpha_DTS_CI, alpha_MK):
    # created by G Oldford Feb 2024
    # purpose: make it easier to get slopes and summary stats dictionary
    #          for a 2D input dataset. Intended for analysis of trends at Nanoose station, model-obs comparisons.

    #    resolution  - float; the precision of the measurement (sig digits) for ties calculated during some stats calcs
    #    remove_nans - True/False;
    #    time_inc    - string; sub-annual time increment for temporal avg ("biweekly", "monthly", "seasonal", "annual")
    #    alpha_CI    - the alpha confidence interval threshold (e.g., 0.05 for 95%)

    # returns:
    #   table                -
    #   summary_stats_all    -
    #   DTS_CI_out           - results of trend test using decorr. time scale
    #                        - list of dict elements; {season: {slope, inter, r_val, p_val, std_err, slope_lcl, slope_ucl}}
    #   mk_3pw_out           - results of '3PW' trend test using Mann-Kendall, see https://doi.org/10.5194/amt-13-6945-2020
    #                        - list of dict elements; {P= max(P_PW, P_TFPW_Y);'ss' (float): statistical significance:
    #                          alpha_MK if the test is ss at the alpha confidence level. Defaults to 95.
    #                          0 if the test is not ss at the alpha_MK confidence level.
    #                          -1 if the test is a TFPW_Y false positive at alpha_MK confidence level
    #                          -2 if the test is a PW false positive at alpha_MK confidence level
    #                          'slope' (float): Sen's slope in units/y
    #                           'ucl' (float): upper confidence level in units/y
    #                           'lcl' (float): lower confidence level in units/y
    #
    #   LR_slopes,
    #   SS_slopes
    #   nan_mask            - return the mask used just in case a 'blind' approach used in subsequent analysis (eg model)


    # slopes from linear regression and alt method theil-sen aka sen's
    LR_slopes = do_lr_seasons(dat_davg, dat_seas, dat_numerictime, dat_timenumeric_seas, seasons_per_year)
    SS_slopes = do_ss_seasons(dat_davg, dat_seas, dat_pytime, dat_pytime_seas, seasons_per_year, resolution, alpha_MK)

    # ==== Detrend ====
    # as in mannkendall lib (mk_white.py)
    slope_method = 'SS' # To do - make "LR" work
    if slope_method == 'SS': slopes = SS_slopes
    elif slope_method == "LR": slopes = LR_slopes
    else: print("error with slope method")
    # note that low p-vals in Het-White and Goldfeld-Quandt tests suggests evidence against homoscedasticity
    summary_stats_all = do_summary_stats(slopes, slope_method, dat_ts, dat_davg, dat_numerictime, dat_seas, dat_timenumeric_seas)

    # ====  Parametric Trend Detection =====
    # use LR or SS with adjusted ESS from DST
    # to-do: make this automatic choice
    type_CI = 'parametric'
    if type_CI == 'parametric':
        DTS_CI_out = dst_trend_sig(summary_stats_all, dat_ts, alpha_DTS_CI, dat_davg, dat_seas, dat_numerictime, dat_timenumeric_seas)

    # ==== Non-Parametric Trend Detection ====
    # - MK w/ Prewhitening (e.g., 3PW)
    type_CI = 'nonparametric'
    if type_CI == 'nonparametric':
        print("MK w/ Sens and 3PW method")
        mk_3pw_out = mk3pw_trend_sig(seasons_per_year, resolution, alpha_MK, dat_davg, dat_pytime, dat_seas, dat_pytime_seas)
    print(mk_3pw_out)

    # write to tables
    # outputs are: LR Slopes, SS Slopes, Summary Stats, MK 3PW Out, DTS CI Out
    headers = list(summary_stats_all[0]['season 1'].keys())
    rows = []
    seas = 1
    for s in summary_stats_all:
        if seas == len(summary_stats_all):
            seas_name = 'all seasons'
            row = s[seas_name].values()
        else:
            seas_name = 'season ' + str(seas)
            row = s[seas_name].values()
        seas += 1
        row = [seas_name] + list(row)
        rows.append(row)
    table = tabulate(rows, headers, tablefmt='grid')

    return (table, summary_stats_all,
            mk_3pw_out, DTS_CI_out,
            LR_slopes, SS_slopes)


# # ==== Paths ====
# meshm_p = '../../data/mesh mask/'
# meshm_f = 'mesh_mask_20210406.nc'
# dat_p = '../climatol_intermediate_files/'
#
# year_min = 1970
# year_max = 2018
# model_or_obs = "Observed" # observed or modelled for labels
# #dat_f = 'Nanoose_obs_anom_temp_1970-2005.nc'
# #dat_f = 'Nanoose_obs_anom_temp_1980-2018.nc'
# dat_f = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
# #dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
# #dat_f = 'RUN203mod_anom_temp_1980-2018.nc'
#
# # ==== Params for Analysis ====
# depth_min = 4.5
# depth_max = 400
# var = 'temperature'
# remove_nans = True
# alpha_DTS_CI = 0.05
# alpha_MK = 0.95
# time_inc = 'monthly' # biweekly, monthly, seasonal, annual
# resolution = 0.0001  # of measurements to use in ties calc in Kendalls Tau and S

#
# (dat_davg, dat_pytime, dat_numerictime, dat_seas,
#  dat_pytime_seas, dat_timenumeric_seas,
#  table, summary_stats_all, mk_3pw_out,
#  DTS_CI_out, LR_slopes, SS_slopes) = do_prelim_trend_analysis(meshm_p, meshm_f, dat_p, dat_f,
#                                                               depth_min, depth_max,
#                                                               var, resolution, remove_nans,
#                                                               alpha_DTS_CI, alpha_MK, time_inc)
#
# # BELOW IS EXPERIMENTAL, AS OF FEB 2024
# # See separate scripts for plots, visuals for manuscript
#
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# from pyemd import EMD, CEEMDAN
# # Apply Empirical Mode Decomposition (EMD)
# emd = EMD()
# IMFs = emd.emd(np.asarray(dat_davg))
# # Plot the original time series and its IMFs
# plt.figure(figsize=(12, 6))
# fig, ax = plt.subplots()
# ax.plot(dat_pytime, dat_davg, label='Original Anomalies')
# for i, imf in enumerate(IMFs):
#     if i > (len(IMFs) - 3):
#         ax.plot(dat_pytime, imf, label=f'IMF {i+1}')
# ax.xaxis.set_major_locator(mdates.YearLocator())
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
#
# plt.xlabel('Time (months)')
# plt.ylabel(str(var) + ' Anomalies')
# plt.title('EMD Decomposition of ' + var + ' Anom. ' + str(year_min) +
#           '-' + str(year_max) + ' ' + str(depth_min) + '-' + str(depth_max))
# plt.legend()
# plt.show()
#
# slope_IMF = IMFs[-1]
# SS_slopes_IMF = do_ss(slope_IMF, dat_pytime, resolution)
# print(SS_slopes_IMF)
#
# ceemdan = CEEMDAN()
# cIMFs = ceemdan.ceemdan(np.asarray(dat_davg))
# slope_cIMF = cIMFs[-1]
# SS_slopes_cIMF = do_ss(slope_cIMF, dat_pytime, resolution)
# print(SS_slopes_cIMF)
#
#
# plt.figure(figsize=(12, 6))
# fig, ax = plt.subplots()
# ax.plot(dat_pytime, dat_davg, label='Original Anomalies')
# for i, imf in enumerate(cIMFs):
#     if i > (len(cIMFs) - 3):
#         if i == len(cIMFs) - 1:
#             label = f'cIMF {i+1}, Slope: ' + str(round(SS_slopes_cIMF[0]['all seasons']['slope'].values()))
#         else:
#             label = f'cIMF {i+1}'
#         ax.plot(dat_pytime, imf, label=label)
#
# plt.xlabel('Time (months)')
# plt.ylabel('Temperature Anomalies')
# ax.xaxis.set_major_locator(mdates.YearLocator())
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
# #
# plt.title(model_or_obs + ' CEEMDAN Decomp of ' + var + ' Anom. ' + str(year_min) +
#           '-' + str(year_max) + ' ' + str(depth_min) + '-' + str(depth_max))
# plt.legend()
# plt.show()




# temp, to get a plot
import matplotlib.pyplot as plt


#dat_f = 'RUN216mod_anom_temp_1980-2018.nc'
# mod = xr.open_dataset(os.path.join(dat_p,dat_f))
# mod = mod[var]
# mod_davg = depth_int(mod, depth_min, depth_max, gdept_0, e3t0)
# mod_davg = mod_davg[~nan_mask]
# mod_davg, _, _ = do_seasonal_avg(mod_davg,time_inc)
#
# fig, ax = plt.subplots()
# ax.fill_between(dat_pytime, dat_davg, mod_davg, color='grey', label='', alpha=0.7, zorder=0)
# #ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
# #ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
# ax.plot(dat_pytime, dat_davg, label='Obs', color='r', linestyle='-', linewidth=1)
# ax.plot(dat_pytime, mod_davg, label='Mod', color='g', linestyle='--', linewidth=1)
# ax.xaxis.set_major_locator(mdates.YearLocator())
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
# ax.set_ylim(-1,1)
# plt.title("HOTSS v1.02 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
# plt.xlabel('Time')
# plt.ylabel('Anomaly (deg C)')
# plt.legend()
# plt.show()
#
# fig, ax = plt.subplots()
# window_size = 5
# mvg_avg_mod = mod_davg.rolling(time_counter=window_size, center=True).mean()
# mvg_avg_obs = dat_davg.rolling(time_counter=window_size, center=True).mean()
# nan_mask = np.isnan(mvg_avg_obs)
# mvg_avg_obs = mvg_avg_obs[~nan_mask]
# mvg_avg_mod = mvg_avg_mod[~nan_mask]
# dat_pytime = dat_pytime[~nan_mask]
# ax.fill_between(dat_pytime, mvg_avg_obs.data, mvg_avg_mod.data, color='grey', label='', alpha=0.7, zorder=0)
# #ax.scatter(dat_pytime, dat_davg, label='Obs', s=1, color='r', marker='o')
# #ax.scatter(dat_pytime, mod_davg, label='Mod', s=1, color='g', marker='x')
# ax.plot(dat_pytime, mvg_avg_obs.data, label='Obs', color='r', linestyle='-', linewidth=1)
# ax.plot(dat_pytime, mvg_avg_mod.data, label='Mod', color='g', linestyle='--', linewidth=1)
# ax.xaxis.set_major_locator(mdates.YearLocator())
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
# ax.set_ylim(-1,1)
# plt.title("HOTSS v1.02 Modelled vs. Observed Temperature Anomalies at Nanoose Stn, " + str(depth_min) + " to " + str(depth_max) + " m")
# plt.xlabel('Time')
# plt.ylabel('Anomaly (deg C)')
# plt.legend()
# plt.show()