# this is a revision of the figure originally done in the python notebook
# higher resolution
# G Oldford - Jan 11, 2025
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import os
import cmocean as cm
from scipy import interpolate


outpath1 = '..//..//data//eval//nanoose_clima//'
mod_run = "RUN216"

mod_salt_bm = xr.open_dataset(outpath1 + 'mod_' + mod_run + '_salinity_1980-2018-bimonthly_timeseries.nc')
mod_temp_bm = xr.open_dataset(outpath1 + 'mod_' + mod_run + '_temperature_1980-2018-bimonthly_timeseries.nc')
obs_salt_clim_modyr = xr.open_dataset(os.path.join(outpath1, 'obs_salinity_1980-2018-climatology.nc'))
obs_temp_clim_modyr = xr.open_dataset(os.path.join(outpath1, 'obs_temperature_1980-2018-climatology.nc'))
mod_salt_clim = xr.open_dataset(outpath1 + 'mod_' + mod_run + '_salinity_1980-2018-climatology.nc')
mod_temp_clim = xr.open_dataset(outpath1 + 'mod_' + mod_run + '_temperature_1980-2018-climatology.nc')

temp_clim = obs_temp_clim_modyr
salt_clim = obs_salt_clim_modyr
obs_temp_clim = obs_temp_clim_modyr
obs_salt_clim = obs_salt_clim_modyr

q_t = (mod_temp_clim - temp_clim).temperature
q_s = (mod_salt_clim - salt_clim).salinity

# plot
fig_dim_1 = 14
fig_dim_2 = 6
grid_rows = 2
grid_cols = 64
grd_hspace = 0.1
grd_wspace = 0.4
fs_x = 8
fs_y = 8
fs1 = 8

fig = plt.figure(figsize=(fig_dim_1, fig_dim_2))
gs = gridspec.GridSpec(grid_rows, grid_cols,
                       hspace=grd_hspace,
                       wspace=grd_wspace)

ax_tclim_obs = plt.subplot(gs[1:2, 5:18])  # bottom left
ax_sclim_obs = plt.subplot(gs[0:1, 5:18], sharex=ax_tclim_obs)  # top left
ax_tclim_mod = plt.subplot(gs[1:2, 18:31], sharey=ax_tclim_obs)  # top middle
ax_sclim_mod = plt.subplot(gs[0:1, 18:31], sharey=ax_sclim_obs, sharex=ax_tclim_mod)  # bottom middle
ax_sdiff = plt.subplot(gs[0:1, 37:50], sharey=ax_sclim_obs, sharex=ax_sclim_mod)  # top right
ax_tdiff = plt.subplot(gs[1:2, 37:50], sharey=ax_tclim_obs, sharex=ax_tclim_mod)  # bottom right

pm_sclim_obs = xr.plot.pcolormesh(obs_salt_clim.salinity, ax=ax_sclim_obs, add_colorbar=False, vmin=19, vmax=32)
pm_tclim_obs = xr.plot.pcolormesh(obs_temp_clim.temperature, ax=ax_tclim_obs, add_colorbar=False, cmap=cm.cm.thermal,
                                  vmin=6, vmax=18)
pm_sclim_mod = xr.plot.pcolormesh(mod_salt_clim.salinity, ax=ax_sclim_mod, add_colorbar=False, vmin=19, vmax=32)
pm_tclim_mod = xr.plot.pcolormesh(mod_temp_clim.temperature, ax=ax_tclim_mod, add_colorbar=False, cmap=cm.cm.thermal,
                                  vmin=6, vmax=18)
pm_sdiff = xr.plot.pcolormesh(q_s, ax=ax_sdiff, add_colorbar=False, cmap=cm.cm.balance)
pm_tdiff = xr.plot.pcolormesh(q_t, ax=ax_tdiff, add_colorbar=False, cmap=cm.cm.balance)

# dummy axes for colorbars
ax_sclim_cb = fig.add_subplot(gs[0:1, 1:2])  # Dummy Axes for colorbars
ax_tclim_cb = fig.add_subplot(gs[1:2, 1:2])
ax_sdiff_cb = fig.add_subplot(gs[0:1, 34:35])
ax_tdiff_cb = fig.add_subplot(gs[1:2, 34:35])

# cbar settings
cb1 = plt.colorbar(pm_sclim_mod, cax=ax_sclim_cb)
cb2 = plt.colorbar(pm_tclim_mod, cax=ax_tclim_cb)
cb3 = plt.colorbar(pm_sdiff, cax=ax_sdiff_cb)
cb4 = plt.colorbar(pm_tdiff, cax=ax_tdiff_cb)
for cb in [cb1, cb2, cb3, cb4]:
    cb.ax.tick_params(labelsize=fs1, left=True, right=False, labelleft=True, labelright=False)
    cb.ax.yaxis.set_label_position("left")

cb1.ax.set_ylabel('Salinity (PSU)', fontdict={'fontsize': 8}, labelpad=0.1)
cb2.ax.set_ylabel('Temperature ($^\circ$C)', fontdict={'fontsize': 8}, labelpad=0.1)
cb3.ax.set_ylabel('Bias (PSU)', fontdict={'fontsize': 8}, labelpad=0.1)
cb4.ax.set_ylabel('Bias ($^\circ$C)', fontdict={'fontsize': 8}, labelpad=0.1)

# first col
for a in [ax_sclim_obs, ax_tclim_obs]:
    plt.sca(a)
    plt.tick_params(axis='y', labelsize=fs_y, pad=1)
    a.set_ylabel('Depth (m)', fontsize=fs1, labelpad=-1)

# second col
for a in [ax_sclim_mod, ax_tclim_mod]:
    plt.sca(a)
    plt.tick_params(axis='y', labelleft=False)
    plt.ylabel('')

# last col
for a in [ax_sdiff, ax_tdiff]:
    plt.sca(a)
    plt.tick_params(axis='y', labelleft=True, labelsize=fs_y)
    plt.ylabel('')
# top row
for a in [ax_sclim_obs, ax_sclim_mod, ax_sdiff]:
    plt.sca(a)
    plt.tick_params(axis='x', labelbottom=False)
    plt.xlabel('')

# bottom row
my_xticks_lab = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
# axs[0].set_xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], my_xticks_lab)
for a in [ax_tclim_obs, ax_tclim_mod, ax_tdiff]:
    plt.sca(a)
    plt.tick_params(axis='x', labelsize=fs_x)
    a.set_xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], my_xticks_lab)
    a.set_xlabel('Month')

# all plots
letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
# let_clrs = ['k', 'k', 'k', 'w', 'w', 'w']
let_clrs = ['k', 'k', 'k', 'k', 'k', 'k']
let_ys = [1.04, 1.04, 1.04, 1.01, 1.01, 1.01]
l = 0
for a in [ax_sclim_obs, ax_sclim_mod, ax_sdiff, ax_tclim_obs, ax_tclim_mod, ax_tdiff]:
    a.set_ylim([1, 400])
    a.set_yscale("log")
    a.invert_yaxis()
    a.set_ylim([400, 2])  # setting to 2 makes the contours better
    letter = letters[l]
    let_clr = let_clrs[l]
    let_y = let_ys[l]
    # a.text(0.05, 0.93, letter, transform=a.transAxes, ha='center', color=let_clr, fontsize=11, zorder=1000)
    a.text(0.05, let_y, letter, transform=a.transAxes, ha='center', color=let_clr, fontsize=10, zorder=1000)
    l += 1

ax_sclim_obs.set_title('Observations')
# ax_sclim_mod.set_title('Model (' + mod_run + ')' )
ax_sclim_mod.set_title('Model')
ax_sdiff.set_title('Model-Observations')

# contours (should combine to one loop)
for s_clim, s_clim_ax in [[obs_salt_clim.salinity, ax_sclim_obs], [mod_salt_clim.salinity, ax_sclim_mod]]:
    X, Y = np.meshgrid(s_clim.timeperiod.values, s_clim.deptht.values)
    levels = [21, 24, 26, 27, 28, 29, 30, 31]
    cont1 = s_clim_ax.contour(X, Y, s_clim, levels=levels, colors='black',
                              linestyles='solid', linewidths=1, alpha=0.7,
                              corner_mask=False
                              )
    s_clim_ax.clabel(cont1, cont1.levels, inline=True, fmt='%.0f', fontsize=8)

for t_clim, t_clim_ax in [[obs_temp_clim.temperature, ax_tclim_obs], [mod_temp_clim.temperature, ax_tclim_mod]]:
    X, Y = np.meshgrid(t_clim.timeperiod.values, t_clim.deptht.values)
    levels = [9, 11, 13, 15, 17, 19]
    cont1 = t_clim_ax.contour(X, Y, t_clim, levels=levels, colors='black',
                              linestyles='solid', linewidths=1, alpha=0.7,
                              corner_mask=False
                              )
    t_clim_ax.clabel(cont1, cont1.levels, inline=True, fmt='%.0f', fontsize=8)

for diff, diff_ax in [[q_s, ax_sdiff],
                      [q_t, ax_tdiff]
                      ]:
    X, Y = np.meshgrid(diff.timeperiod.values, diff.deptht.values)
    #     levels=[9,11,13,15,17,19]
    cont1 = diff_ax.contour(X, Y, diff, colors='black',
                            linestyles='solid', linewidths=1, alpha=0.7,
                            corner_mask=False
                            )
    diff_ax.clabel(cont1, cont1.levels, inline=True, fmt='%.0f', fontsize=8)

#plt.savefig('../../figs/Fig08_rev011125.png', dpi=600)
#plt.savefig('../../figs/Fig08_rev011125.pdf', dpi=600)