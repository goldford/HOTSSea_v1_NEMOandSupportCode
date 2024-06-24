# Created Feb 9 2024 by G Oldford
# Purpose: seasonal and nonseasonal trend estimation over all depth levs
# python 3.11.5 - set the interpeter to this (anaconda3 default used here)
# statsmodels v.0.14, netCDF4 v1.6.2, mannkendall 1.1.1, matplotlib v3.6.3

from GO_tools import get_dat, do_ss_seasons
from statsmodels.tsa.stattools import acf
import numpy as np
import mannkendall as mk
import os
import matplotlib.pyplot as plt
import netCDF4 as nc

# ==== Paths ====
meshm_p = '..//..//data//mesh mask//'
meshm_f = 'mesh_mask_20210406.nc'
dat_p = '..//..//data//eval//nanoose_clima//'

# to-do: this does not currently do much - just selects files below
year_min = 1980
year_max = 2018

#dat_f = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
#dat_f = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
dat_f_obs = 'Nanoose_obs_anom_temp_' + str(year_min) + '-' + str(year_max) + '.nc'
dat_f_obs_full = 'Nanoose_obs_anom_temp_' + str(1970) + '-' + str(year_max) + '.nc'
dat_f_mod1 = 'RUN203mod_anom_temp_1980-2018.nc'
dat_f_mod2 = 'RUN216mod_anom_temp_1980-2018.nc'

# ==== Params for Analysis ====
var = 'temperature'
remove_nans = True
use_abs = False # make sure false if you want trends in variable (not variability)!!
alpha_DTS_CI = 0.05
alpha_MK = 95
time_inc = 'seasonal' # biweekly, monthly, seasonal, annual
resolution = 0.0001  # of measurements to use in ties calc in Kendalls Tau and S
if time_inc == 'annual': seasons_per_year = 1
elif time_inc == 'seasonal': seasons_per_year = 4
elif time_inc == 'monthly': seasons_per_year = 12
elif time_inc == 'biweekly': seasons_per_year = 24

# trend for all depths - experimental - meant for obs_full
yr_st = 0
yr_en = 0 # change to 0 not to do experiment
(d_o, d_pt_o, d_tn_o,
 d_se_o, d_pt_se_o, d_tn_se_o,
 d_ts, time_dim, se_year) = get_dat(meshm_p, meshm_f, dat_p, dat_f_obs,
                                    var, time_inc, use_abs, dep_int=False, yr_st=yr_st, yr_en=yr_en)

# slopes
SS_slopes = do_ss_seasons(d_o, d_se_o, d_pt_o, d_pt_se_o, seasons_per_year, resolution, alpha_MK)

# test for autocorr in slope residuals
# autocorrelation
alpha_acf=0.05
cnt_mk3 = 0
mk3_on = True
for season in range(0, seasons_per_year):
    d_np = np.asarray(d_se_o[season])
    d_tn = d_tn_se_o[season] # time numeric
    d_pt = d_pt_se_o[season] # python datetime
    for d_idx in range(d_np.shape[0]):
        depth_data = d_np[d_idx, :]
        # get the slope and intercept from the dictionary
        # awkward... can streamline somehow? - GO 20240209
        i = 0
        for sr in SS_slopes:
            for s in sr:
                snm = "season " + str(season+1)
                if s == snm:
                    for key, val in sr.items():
                        d_key = 'depth idx' + str(d_idx)
                        if d_key in val:
                            slope = val[d_key]['slope']
                            inter = val[d_key]['inter']

                            residuals = depth_data - (slope * (d_tn / 365) + inter)
                            nlags = len(residuals - 1)
                            acf_values = acf(residuals.values, adjusted=True, qstat=True, fft=False, alpha=alpha_acf,
                                             nlags=nlags,
                                             missing="conservative")
                            AR1_coef = acf_values[0][1]
                            AR1_confint = acf_values[1]
                            AR1_qstat = acf_values[2]
                            AR1_pval = acf_values[3]  # note low vals means autocorr
                            SS_slopes[i][s][d_key]['acf-ar1-coeff'] = AR1_coef
                            SS_slopes[i][s][d_key]['acf-qstat'] = AR1_qstat
                            SS_slopes[i][s][d_key]['acf-pvals'] = AR1_pval
                            SS_slopes[i][s][d_key]['acf-pval_1'] = AR1_pval[0]
                            SS_slopes[i][s][d_key]['acf-all-nlag-coeffs'] = acf_values

                            # if pval > alpha then recompute slope, CI, etc using mk 3pw
                            if mk3_on:
                                if AR1_pval[0] <= alpha_acf:
                                    out = mk.mk_temp_aggr(d_pt, depth_data, resolution)
                                    SS_slopes[i][s][d_key]['slope'] = out[0]['slope']
                                    SS_slopes[i][s][d_key]['ucl'] = out[0]['ucl']
                                    SS_slopes[i][s][d_key]['lcl'] = out[0]['lcl']
                                    SS_slopes[i][s][d_key]['p'] = out[0]['p']
                                    inter_yr = np.nanmedian(depth_data) - np.median(
                                        np.arange(len(depth_data))[~np.isnan(depth_data.flatten())]) * out[0]['slope']
                                    SS_slopes[i][s][d_key]['inter'] = inter_yr
                                    SS_slopes[i][s][d_key]['ss'] = out[0]['ss']
                                    cnt_mk3 += 1
            i += 1

# repeat but for 'all seasons' together
d = d_o
d_tn = d_tn_o
d_np = np.asarray(d)
d_pt = d_pt_o
for d_idx in range(d_np.shape[0]):
    depth_data = d_np[d_idx, :]
    # get the slope and intercept from the dictionary
    # awkward... can streamline somehow? - GO 20240209
    i = 0
    for sr in SS_slopes:
        for s in sr:
            snm = "all seasons"
            if s == snm:
                for key, val in sr.items():
                    d_key = 'depth idx' + str(d_idx)
                    if d_key in val:
                        slope = val[d_key]['slope']
                        inter = val[d_key]['inter']
                        residuals = depth_data - (slope * (d_tn / 365) + inter)
                        nlags = len(residuals - 1)
                        acf_values = acf(residuals.values, adjusted=True, qstat=True, fft=False, alpha=alpha_acf,
                                         nlags=nlags,
                                         missing="conservative")
                        AR1_coef = acf_values[0][1]
                        AR1_confint = acf_values[1]
                        AR1_qstat = acf_values[2]
                        AR1_pval = acf_values[3]  # note low vals means autocorr
                        SS_slopes[i][s][d_key]['acf-ar1-coeff'] = AR1_coef
                        SS_slopes[i][s][d_key]['acf-qstat'] = AR1_qstat
                        SS_slopes[i][s][d_key]['acf-pvals'] = AR1_pval
                        SS_slopes[i][s][d_key]['acf-pval_1'] = AR1_pval[0]
                        SS_slopes[i][s][d_key]['acf-all-nlag-coeffs'] = acf_values

                        # if pval > alpha then recompute slope, CI, etc using mk 3pw
                        if mk3_on:
                            if AR1_pval[0] <= alpha_acf:
                                out = mk.mk_temp_aggr(d_pt, depth_data, resolution)
                                SS_slopes[i][s][d_key]['slope'] = out[0]['slope']
                                SS_slopes[i][s][d_key]['ucl'] = out[0]['ucl']
                                SS_slopes[i][s][d_key]['lcl'] = out[0]['lcl']
                                SS_slopes[i][s][d_key]['p'] = out[0]['p']
                                inter_yr = np.nanmedian(depth_data) - np.median(
                                    np.arange(len(depth_data))[~np.isnan(depth_data.flatten())]) * out[0]['slope']
                                SS_slopes[i][s][d_key]['inter'] = inter_yr
                                SS_slopes[i][s][d_key]['ss'] = out[0]['ss']
                                cnt_mk3 += 1
        i += 1

print('number of re-analyses using mk3 due to autocorre:', cnt_mk3)
print('out of a total of ', len(SS_slopes) * d_np.shape[0])
SS_slopes_o = SS_slopes


############################################################
############################################################
# REPEAT for model
# should revise to reduce this repetition
#yr_st=0 # set to zero for default
#yr_en=0
#(d_m, d_pt_m, d_tn_m,
# d_se_m, d_pt_se_m, d_tn_se_m, _, _, _) = get_dat(meshm_p, meshm_f, dat_p, dat_f_obs, var,
#                                                     time_inc, use_abs, dep_int=False, yr_st=yr_st, yr_en=yr_en)
yr_st=0
yr_en=0
(d_m, d_pt_m, d_tn_m,
 d_se_m, d_pt_se_m, d_tn_se_m, _, _, _) = get_dat(meshm_p, meshm_f, dat_p, dat_f_mod2, var,
                                                      time_inc, use_abs, dep_int=False, yr_st=yr_st, yr_en=yr_en)

# slopes
SS_slopes = do_ss_seasons(d_m, d_se_m, d_pt_m, d_pt_se_m, seasons_per_year, resolution, alpha_MK)
# test for autocorr in slope residuals
# autocorrelation
alpha_acf=0.05
cnt_mk3 = 0
mk3_on = True
for season in range(0, seasons_per_year):
    d_np = np.asarray(d_se_m[season])
    d_tn = d_tn_se_m[season] # time numeric
    d_pt = d_pt_se_m[season] # python datetime
    for d_idx in range(d_np.shape[0]):
        depth_data = d_np[d_idx, :]
        # get the slope and intercept from the dictionary
        # awkward... can streamline somehow? - GO 20240209
        i = 0
        for sr in SS_slopes:
            for s in sr:
                snm = "season " + str(season+1)
                if s == snm:
                    for key, val in sr.items():
                        d_key = 'depth idx' + str(d_idx)
                        if d_key in val:
                            slope = val[d_key]['slope']
                            inter = val[d_key]['inter']

                            residuals = depth_data - (slope * (d_tn / 365) + inter)
                            nlags = len(residuals - 1)
                            acf_values = acf(residuals.values, adjusted=True, qstat=True, fft=False, alpha=alpha_acf,
                                             nlags=nlags,
                                             missing="conservative")
                            AR1_coef = acf_values[0][1]
                            AR1_confint = acf_values[1]
                            AR1_qstat = acf_values[2]
                            AR1_pval = acf_values[3]  # note low vals means autocorr
                            SS_slopes[i][s][d_key]['acf-ar1-coeff'] = AR1_coef
                            SS_slopes[i][s][d_key]['acf-qstat'] = AR1_qstat
                            SS_slopes[i][s][d_key]['acf-pvals'] = AR1_pval
                            SS_slopes[i][s][d_key]['acf-pval_1'] = AR1_pval[0]
                            SS_slopes[i][s][d_key]['acf-all-nlag-coeffs'] = acf_values

                            # if pval > alpha then recompute slope, CI, etc using mk 3pw
                            if mk3_on:
                                if AR1_pval[0] <= alpha_acf:
                                    out = mk.mk_temp_aggr(d_pt, depth_data, resolution)
                                    SS_slopes[i][s][d_key]['slope'] = out[0]['slope']
                                    SS_slopes[i][s][d_key]['ucl'] = out[0]['ucl']
                                    SS_slopes[i][s][d_key]['lcl'] = out[0]['lcl']
                                    SS_slopes[i][s][d_key]['p'] = out[0]['p']
                                    inter_yr = np.nanmedian(depth_data) - np.median(
                                        np.arange(len(depth_data))[~np.isnan(depth_data.flatten())]) * out[0]['slope']
                                    SS_slopes[i][s][d_key]['inter'] = inter_yr
                                    SS_slopes[i][s][d_key]['ss'] = out[0]['ss']
                                    cnt_mk3 += 1
            i += 1

# repeat but for 'all seasons' together
d = d_m
d_tn = d_tn_m
d_np = np.asarray(d)
d_pt = d_pt_m
for d_idx in range(d_np.shape[0]):
    depth_data = d_np[d_idx, :]
    # get the slope and intercept from the dictionary
    # awkward... can streamline somehow? - GO 20240209
    i = 0
    for sr in SS_slopes:
        for s in sr:
            snm = "all seasons"
            if s == snm:
                for key, val in sr.items():
                    d_key = 'depth idx' + str(d_idx)
                    if d_key in val:
                        slope = val[d_key]['slope']
                        inter = val[d_key]['inter']
                        residuals = depth_data - (slope * (d_tn / 365) + inter)
                        nlags = len(residuals - 1)
                        acf_values = acf(residuals.values, adjusted=True, qstat=True, fft=False, alpha=alpha_acf,
                                         nlags=nlags,
                                         missing="conservative")
                        AR1_coef = acf_values[0][1]
                        AR1_confint = acf_values[1]
                        AR1_qstat = acf_values[2]
                        AR1_pval = acf_values[3]  # note low vals means autocorr
                        SS_slopes[i][s][d_key]['acf-ar1-coeff'] = AR1_coef
                        SS_slopes[i][s][d_key]['acf-qstat'] = AR1_qstat
                        SS_slopes[i][s][d_key]['acf-pvals'] = AR1_pval
                        SS_slopes[i][s][d_key]['acf-pval_1'] = AR1_pval[0]
                        SS_slopes[i][s][d_key]['acf-all-nlag-coeffs'] = acf_values

                        # if pval > alpha then recompute slope, CI, etc using mk 3pw
                        if mk3_on:
                            if AR1_pval[0] <= alpha_acf:
                                out = mk.mk_temp_aggr(d_pt, depth_data, resolution)
                                SS_slopes[i][s][d_key]['slope'] = out[0]['slope']
                                SS_slopes[i][s][d_key]['ucl'] = out[0]['ucl']
                                SS_slopes[i][s][d_key]['lcl'] = out[0]['lcl']
                                SS_slopes[i][s][d_key]['p'] = out[0]['p']
                                inter_yr = np.nanmedian(depth_data) - np.median(
                                    np.arange(len(depth_data))[~np.isnan(depth_data.flatten())]) * out[0]['slope']
                                SS_slopes[i][s][d_key]['inter'] = inter_yr
                                SS_slopes[i][s][d_key]['ss'] = out[0]['ss']
                                cnt_mk3 += 1
        i += 1

print('number of re-analyses using mk3 due to autocorre:', cnt_mk3)
print('out of a total of ', len(SS_slopes) * d_np.shape[0])
# if the autocorr p value

SS_slopes_m = SS_slopes


# plot over depths
season_names = ['season 1', 'season 2', 'season 3', 'season 4', 'inter-seasonal'] # 'all seasons' is not good to use
season_labels = ['Winter', 'Spring', 'Summer', 'Fall', 'Annual']
letters = ['(a)','(b)','(c)','(d)','(e)']
fig, axs = plt.subplots(1, len(season_names), figsize=(7, 3.5), sharey=True)
log_on = True

for i, season_name in enumerate(season_names):

    # Omitting the annualised trend
    # if i == len(season_names) - 1:
    #     continue
    season_label = season_labels[i]

    #############################################
    #############################################
    # Observations

    depidx = []
    depths = []
    trends = []
    ucls = []
    lcls = []
    ps = []

    # Extract trend, ucl, and lcl for each depth for the specified season
    for sr in SS_slopes_o:
        for s in sr:
            snm = season_name
            if s == snm:
                for key, val in sr.items():
                    for key2, val2 in val.items():
                        depidx.append(val2['depth idx'])
                        trends.append(val2['slope']*10)
                        ucls.append(val2['ucl']*10)
                        lcls.append(val2['lcl']*10)
                        ps.append(val2['p'])
    # get real depths, not the idx's
    depidx = np.array(depidx)
    row_j = 175; col_i = 75
    with nc.Dataset(os.path.join(meshm_p, meshm_f), 'r') as mesh:
        gdept_0 = mesh.variables['gdept_0'][:]
        gdept_0 = gdept_0[0, :, row_j, col_i]
    depths = []
    for didx in depidx: depths.append(gdept_0[didx])
    trends = np.array(trends)
    ucls = np.array(ucls)
    lcls = np.array(lcls)
    ps = np.array(ps)

    # Plotting
    lbl='Obs.'
    clr = 'k'
    #lbl = '1980 - 2018'
    #clr = 'k'
    if log_on:
        trends = trends[:]
        depths = depths[:]
        ucls = ucls[:]
        lcls = lcls[:]
    axs[i].plot(trends, depths, label=lbl, color=clr, linewidth=1.2)
    #axs[i].plot(ucls, depths, linestyle='--', label='', color=clr, linewidth=1)
    #axs[i].plot(lcls, depths, linestyle='--', label='', color=clr, linewidth=1)
    axs[i].fill_betweenx(depths, ucls, lcls, color=clr, alpha=0.2, linewidth=0)

    #########################################################
    #########################################################
    # Repeat for Model
    depidx = []
    depths = []
    trends = []
    ucls = []
    lcls = []
    ps = []

    # Extract trend, ucl, and lcl for each depth for the specified season
    for sr in SS_slopes_m:
        for s in sr:
            snm = season_name
            if s == snm:
                for key, val in sr.items():
                    for key2, val2 in val.items():
                        depidx.append(val2['depth idx'])
                        trends.append(val2['slope'] * 10)
                        ucls.append(val2['ucl'] * 10)
                        lcls.append(val2['lcl'] * 10)
                        ps.append(val2['p'])
    # get real depths, not the idx's
    depidx = np.array(depidx)
    row_j = 175; col_i = 75
    with nc.Dataset(os.path.join(meshm_p, meshm_f), 'r') as mesh:
        gdept_0 = mesh.variables['gdept_0'][:]
        gdept_0 = gdept_0[0, :, row_j, col_i]
    depths = []
    for didx in depidx: depths.append(gdept_0[didx])
    trends = np.array(trends)
    ucls = np.array(ucls)
    lcls = np.array(lcls)
    ps = np.array(ps)

    # Plotting
    clr = 'r'
    lbl = 'HOTSSea'
    # clr = 'grey'
    # lbl = '1970 - 1990'
    if log_on:
        trends = trends[:]
        depths = depths[:]
        ucls = ucls[:]
        lcls = lcls[:]
    axs[i].plot(trends, depths, label=lbl, color=clr, linewidth=1.2)
    #axs[i].plot(ucls, depths, linestyle='--', label='', color=clr, linewidth=1)
    #axs[i].plot(lcls, depths, linestyle='--', label='', color=clr, linewidth=1)
    axs[i].fill_betweenx(depths, ucls, lcls, color=clr, alpha=0.2, linewidth=0)

    # if 1975 - 1985
    #xl_min = -1
    #xl_max = 2.5
    # if 1970 - 2018
    #xl_min = -0.5
    #xl_max = 1.25
    # if mod obs
    xl_min = -0.25
    xl_max = 0.6
    axs[i].set_xlim(xl_min, xl_max)
    #axs[i].set_ylim(1, max(depths))

    #axs[i].set_xticks(np.arange(xl_min, xl_max, 0.5)) # if 1975 - 1985
    axs[i].set_xticks(np.arange(xl_min, 0.75, 0.25))

    axs[i].set_title(f'{season_label}', fontsize=8.5)
    axs[i].tick_params(axis='both', which='major', labelsize=6)

    letter = letters[i]
    axs[i].text(0.1, 1.03, letter, transform=axs[i].transAxes, ha='center', color='k', fontsize=9)

    if log_on:
        axs[i].set_yscale("log")

    if i == 0:
        axs[i].set_ylabel('Depth (m)', fontsize=8)
    if i == 2:
        axs[i].set_xlabel('Temp. Anom. Trend (\u00B0C) / decade)', fontsize=8)
    if not log_on and i == 2:
        axs[i].legend(fontsize=6, loc='lower right')
    elif i == 2:
        axs[i].legend(fontsize=6, loc='lower right')

    axs[i].invert_yaxis()  # Invert y-axis to have depths increasing downwards

    axs[i].grid(color='lightgrey')
    axs[i].axvline(x=0, color='black', linestyle='--', linewidth='0.5')


plt.tight_layout()

# plt.savefig('trend_alldep_nanoose_modobs_wletters.png', dpi=300)
plt.savefig('..//..//figs//Fig10_' +
            str(year_min) + '-' + str(year_max) + '.png',
            dpi=300)
# plt.savefig('..//..//figs//Fig10_' +
#             str(year_min) + '-' + str(year_max) + '.eps',
#             dpi=300)
plt.savefig('..//..//figs//Fig10_' +
            str(year_min) + '-' + str(year_max) + '.pdf',
            dpi=300)
plt.show()


print("done")