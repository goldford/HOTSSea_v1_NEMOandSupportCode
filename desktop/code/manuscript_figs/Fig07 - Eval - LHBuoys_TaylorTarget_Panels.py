# Created by G Oldford Aug 2023
# last edit March 2024

# Purpose: Taylor and Target plots of SST and SSS from buoys and lighthouses
# Inputs: buoy and LH analysis .pickle file exported from the 'pyap' package
# Process: loop over the model runs we wish to visualize,
#          optional: use arrows to visualize the effect of using a moving average
#          (presumably to improve fits due to obs bias or error and / or model issues
# Output: panel plot using matplotlib gridspec of results. relies on skillmetrics package

# Notes / To-dos:
#   20230821 - two loops are inefficient (on for each of taylor and target).
#              b/c of opening and re-opening and processing pickles repeatedly
#   20240320 - using skillmetrics 1.1.8 - newer does not work - but this requires python 3.9
#            - ideally don't use skillmetrics in future
#            - numpy 1.21.5, py 3.9.13, matplotlib 3.5.2
#            - skillmetrics used here was edited - get_target_diagram_options,
#            overlay_target_diagram_circles, plot_pattern_diagram_markers,
#            overlay_taylor_diagram_circles, plot_taylor_axes
#            - set to py 3.9 venv before running - projects -> settings etc.


import numpy as np
# attempt to fix issue with skillmetrics library 1.1.8 and py 3.11 - 2024-03 GO


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import skill_metrics as sm
import pickle
import os

from utils_GO import get_taylor_stats_SST, get_taylor_stats_LH, get_target_stats_SST, get_target_stats_LH

## streamline the above
modelruns_info = {  # 'SalishSea1500-RUN203': {'path': 'D:/temp_nemo/RUN203_PLOTS_SCORES/',
    #                                             'colour': 'k',
    #                                             'shortcode': 'RUN203', 'experiment': False},
    #                        'SalishSea1500-RUN216': {'path': 'D:/temp_nemo/RUN216/',
    #                                              'colour': 'k',
    #                                             'shortcode': 'RUN216', 'experiment': False},
    # 'SalishSea1500-RUN216-altSST': {'path': 'D:/temp_nemo/RUN216_altSST/',
    #                                 'colour': 'r', 'colour_face': 'r', 'colour_edge': 'r',
    #                                 'shortcode': 'RUN216', 'experiment': False}  # ,
    # path to old PC above
    # alt SST = z = 1 rather than z = 0
    'SalishSea1500-RUN216-altSST': {'path': 'C:/Users/Greig/Sync/1. UBC-Laptop/5. PhD Dissertation/2_NEMO_SS1500/pyap_out/RUN216_altSST/',
                                    'colour': 'r', 'colour_face': 'r', 'colour_edge': 'r',
                                    'shortcode': 'RUN216', 'experiment': False}  # ,
    #                      'SalishSea500-201905': {'path': 'D:/temp_nemo/SS500/',
    #                                            'colour': 'g',
    #                                            'shortcode': '201905', 'experiment': True}
}


arrow_mutatescale = 11
arrow_facecolor = '#f56f6f'
arrow_edgecolor = '#f56f6f'
arrow_lw = 1

# marker and plotting options
msize1 = 6
alpha = 0.8
OPTIONS_M_SST = {
    # SST / wave buoys (others not used are ptaw1, ptww1,tcnw1,neaw1,46087,46088)
    '46131': {'marker': "o", 'size': msize1, 'label': '46131', 'mf_colour': 'w', 'me_colour': '#000000'},
    '46146': {'marker': "o", 'size': msize1, 'label': '46146', 'mf_colour': 'k', 'me_colour': '#000000'},
    '46134': {'marker': "o", 'size': msize1, 'label': '46134', 'mf_colour': '#814f8d', 'me_colour': '#000000'}
}
# OPTIONS_M_SST = {
#     # SST / wave buoys (others not used are ptaw1, ptww1,tcnw1,neaw1,46087,46088)
#     '46131': {'marker': "$1$", 'size': msize, 'label': '46131'},
#     '46146': {'marker': "$3$", 'size': msize, 'label': '46146'},
#     '46134': {'marker': "$2$", 'size': msize, 'label': '46134'}
# }
msize2=9
OPTIONS_M_LH = {'active_pass_LH.nc': {'marker': 'P', 'size': msize2, 'label': 'Active Pass',
                                      'mf_colour': 'w', 'me_colour': '#000000'},
                'cape_mudge_LH.nc': {'marker': 'P', 'size': msize2, 'label': 'Cape Mudge',
                                     'mf_colour': '#33a02c', 'me_colour': '#000000'},
                'chrome_island_LH.nc': {'marker': 'P', 'size': msize2, 'label': 'Chrome Is.',
                                        'mf_colour': '#000000', 'me_colour': '#000000'},
                'entrance_island_LH.nc': {'marker': 'X', 'size': msize2, 'label': 'Entrance Is.',
                                          'mf_colour': 'w', 'me_colour': '#000000'},
                'race_rocks_LH.nc': {'marker': 'X', 'size': msize2, 'label': 'Race Rocks',
                                     'mf_colour': '#fdb142', 'me_colour': '#000000'},
                'sheringham_point_LH.nc': {'marker': 'X', 'size': msize2, 'label': 'Sheringham Pt.',
                                           'mf_colour': '#000000', 'me_colour': '#000000'},
                'sisters_islets_LH.nc': {'marker': 'X', 'size': msize2, 'label': 'Sisters Islets',
                                         'mf_colour': '#409c5a', 'me_colour': '#000000'}
                }
# OPTIONS_M_LH = {'active_pass_LH.nc': {'marker': r'$\mathdefault{A}$', 'size': msize, 'label': 'Active Pass'},
#                 'cape_mudge_LH.nc': {'marker': r'$\mathdefault{B}$', 'size': msize, 'label': 'Cape Mudge'},
#                 'chrome_island_LH.nc': {'marker': r'$\mathdefault{C}$', 'size': msize, 'label': 'Chrome Is.'},
#                 'entrance_island_LH.nc': {'marker': r'$\mathdefault{D}$', 'size': msize, 'label': 'Entrance Is.'},
#                 'race_rocks_LH.nc': {'marker': r'$\mathdefault{E}$', 'size': msize, 'label': 'Race Rocks'},
#                 'sheringham_point_LH.nc': {'marker': r'$\mathdefault{F}$', 'size': msize, 'label': 'Sheringham Pt.'},
#                 'sisters_islets_LH.nc': {'marker': r'$\mathdefault{G}$', 'size': msize, 'label': 'Sisters Islets'}
#                 }
lh_list = []
for lh in OPTIONS_M_LH.keys():
    lh_list.append(lh)
buoy_list = ['46131', '46134', '46146']  # excluding many of these
mv_avg_list = {1: {'color': 'r'}, 30: {'color': 'k'}}
augment_rmsd = True
use_nanmean = True
use_arrows = True

fig = plt.figure(figsize=(8, 8))

# it is a 2x2 plot, gs gives space for extra spacing col and row for titles, legend
gs = gridspec.GridSpec(13, 12, figure=fig)

# plot the taylor diagrams in first col
var_codes = ["T", "S"]
for var_ts in var_codes:
    # determine which subplot to put plot
    if var_ts == "T":
        ax = fig.add_subplot(gs[1:6, 0:5])
        # rincSTD = 0.25  # increment, not working
        tickSTD = [0.0, 0.5, 1.0, 1.5]
        axismax = 1.5
        tickRMS = [0.5, 1.0, 1.5]
        ax = plt.gca()
        ax.text(-0.08, 1.75, "(a)", ha="left", va="top", fontsize=12)
        ax.text(0.75, 1.75, "Temperature", ha="center", va="top", fontsize=14)
        titleRMS = 'on'

    elif var_ts == 'S':
        ax = fig.add_subplot(gs[1:6, 6:11])
        # rincSTD = 0.25  # increment
        tickSTD = [0.0, 0.5, 1.0, 1.5, 2.0]
        tickRMS = [0.5, 1.0, 1.5, 2.0]
        # ticksCOR = [0.2,0.4,0.6,0.]
        axismax = 2.5
        ax = plt.gca()
        ax.text(-0.08, 2.9, "(b)", ha="left", va="top", fontsize=12)
        ax.text(1.25, 2.9, "Salinity", ha="center", va="top", fontsize=14)
        # ax.text(0.75, -0.3, "Salinity", ha="center", va="top", fontsize=14)
        titleRMS = 'off'

    # the sm.taylor_diagram code uses the current active axis
    sm.taylor_diagram(1, 0, 1, numberPanels=1, axismax=axismax,
                      styleCOR='-', colCOR='#848484',widthcor=0.8,
                      styleSTD='-', colSTD='#848484', tickSTD=tickSTD,  widthstd=0.8,# rincSTD=rincSTD,
                      styleRMS='--', colRMS='#848484', titleRMS=titleRMS, labelRMS='NCRMSE', widthRMS=0.8,
                      rmsLabelFormat='0:.1f', tickRMS=tickRMS, titlermsdangle=120, # determins where on circle label goes
                      markerOBS='*', colOBS='k', markerSize=msize1,
                      overlay='off')

    for run in modelruns_info.keys():
        # hack
        if run == 'SalishSea1500-RUN216-altSST':
            run_sname = 'SalishSea1500-RUN216'
        else:
            run_sname = run

        if var_ts == "T":
            print("Taylor SST", run)
        elif var_ts == "S":
            print("Taylor SSS", run)
        else:
            print("No variable chosen")
            break

        path1 = modelruns_info[run]['path']
        shortcode = modelruns_info[run]['shortcode']

        # plot LH data
        print("plotting LH data")
        LH_pic_data = 'LH_class4_' + shortcode + '_hindcast.pickle'
        with open(os.path.join(path1, LH_pic_data), 'rb') as f:
            LH_scores = pickle.load(f)

        # get stats and plot each lh
        for lh in LH_scores['obs'].keys():
            if lh not in lh_list:
                continue
            print(lh)

            lh_befores = ()
            lh_afters = ()

            plt_n = 0
            for mv_avg in mv_avg_list.keys():
                print("plotting moving average " + str(mv_avg) + ' days')
                color_mvg_avg = mv_avg_list[mv_avg]['color']

                stats_pack = get_taylor_stats_LH(lh, mv_avg,
                                                 LH_scores, run_sname,
                                                 use_nanmean, var_ts)

                # plot a marker
                try:
                    marker_dict = OPTIONS_M_LH[lh]
                    mf_colour = marker_dict['mf_colour']
                    me_colour = marker_dict['me_colour']
                except:
                    print("No marker found for " + lh + ". ")
                    continue
                # if pointing to improvement after mvg avg plot arrows
                if (use_arrows and plt_n > 0):
                    stats_end = stats_pack

                    polar_CORS = np.arccos(np.asarray((stats_start[6])))
                    polar_CORS_after = np.arccos(np.asarray((stats_end[6])))
                    #
                    theta1 = polar_CORS
                    theta2 = polar_CORS_after
                    r1 = np.asarray((stats_start[3]))
                    r2 = np.asarray((stats_end[3]))

                    # Convert polar to Cartesian coordinates
                    x1 = r1 * np.cos(theta1)
                    y1 = r1 * np.sin(theta1)

                    x2 = r2 * np.cos(theta2)
                    y2 = r2 * np.sin(theta2)

                    ax = plt.gca()
                    if mf_colour == 'w':
                        mf_colour = 'k'
                    ax.annotate('', xy=(x2, y2),
                                xytext=(x1, y1),
                                arrowprops=dict(
                                    arrowstyle='-|>',
                                    #head_width=0.7,
                                    #shrink_factor=0.2,
                                    mutation_scale=arrow_mutatescale,
                                    facecolor=mf_colour,
                                    edgecolor=mf_colour,
                                    lw=arrow_lw,
                                    ls='--',
                                    zorder=1
                                ),
                                # arrowprops=dict(
                                #     facecolor='k',
                                #     edgecolor='k',
                                #     lw=0.1,
                                #     ls='--',
                                #     headwidth=3,
                                #     headlength=6
                                # ),
                                fontsize=10,
                                ha='center',
                                va='center'
                                )
                else:
                    stats_start = stats_pack
                    # print("stats for lh " + lh + " "  + str(stats_pack))
                    sm.taylor_diagram(
                        np.asarray((1, stats_pack[3])),  # mod_norm_stdev
                        np.asarray((1, stats_pack[5])),  # ncrmsd
                        np.asarray((1, stats_pack[6])),  # ccoef
                        # note change required to plot_pattern_diagram_markers.py in sm to get
                        # symbols that are letters or numbers to work
                        markersymbol=marker_dict['marker'],
                        # markercolor=modelruns_info[run]['colour'],
                        # markercolor=color_mvg_avg,
                        markercolor=mf_colour,
                        markersize=marker_dict['size'],
                        alpha=alpha,
                        overlay='on'
                    )
                plt_n += 1
        LH_scores = None

        # Plot SST from buoys
        if var_ts == 'T':
            print("plotting buoys..")
            SSTbuoy_pic_data = "SST_class4_" + run + ".pickle"
            with open(os.path.join(path1, SSTbuoy_pic_data), 'rb') as f:
                SSTbuoy_scores = pickle.load(f)

            for bu in SSTbuoy_scores['obs'].keys():

                # hack
                if bu == "C46134_2001-02_2016-12.nc":
                    buoy_bettername = '46134'
                elif bu == "C46131_1992-10_2021-01.nc":
                    buoy_bettername = '46131'
                elif bu == "C46146_1992-03_2022-06.nc":
                    buoy_bettername = '46146'
                elif bu == "C46182_1989-09_1991-11.nc":
                    buoy_bettername = '46182'
                else:
                    buoy_bettername = bu

                if (buoy_bettername not in buoy_list):
                    continue
                print(buoy_bettername)

                plt_n = 0
                for mv_avg in mv_avg_list.keys():
                    print("plotting moving average " + str(mv_avg) + ' days')
                    color_mvg_avg = mv_avg_list[mv_avg]['color']

                    # returns buoy_bettername, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss
                    stats_pack = get_taylor_stats_SST(bu, mv_avg, SSTbuoy_scores,
                                                      run_sname, use_nanmean)
                    # plot one marker
                    try:
                        marker_dict = OPTIONS_M_SST[buoy_bettername]
                        mf_colour = marker_dict['mf_colour']
                        me_colour = marker_dict['me_colour']
                    except:
                        print("No marker found for " + buoy_bettername + ". ")
                        continue

                    # if pointing to improvement after mvg avg plot arrows
                    print(plt_n)
                    if (use_arrows and plt_n > 0):
                        stats_end = stats_pack

                        polar_CORS = np.arccos(np.asarray((stats_start[6])))
                        polar_CORS_after = np.arccos(np.asarray((stats_end[6])))
                        #
                        theta1 = polar_CORS
                        theta2 = polar_CORS_after
                        r1 = np.asarray((stats_start[3]))
                        r2 = np.asarray((stats_end[3]))

                        # Convert polar to Cartesian coordinates
                        x1 = r1 * np.cos(theta1)
                        y1 = r1 * np.sin(theta1)

                        x2 = r2 * np.cos(theta2)
                        y2 = r2 * np.sin(theta2)

                        ax = plt.gca()
                        if mf_colour == 'w':
                            mf_colour = '#585858'
                        ax.annotate('', xy=(x2, y2),
                                    xytext=(x1, y1),
                                    arrowprops=dict(
                                        arrowstyle='-|>',
                                        # head_width=0.7,
                                        # shrink_factor=0.2,
                                        mutation_scale=arrow_mutatescale,
                                        facecolor=mf_colour,
                                        edgecolor=mf_colour,
                                        lw=arrow_lw,
                                        ls='--',
                                        zorder=1
                                    ),
                                    # arrowprops=dict(
                                    #     facecolor='k',
                                    #     edgecolor='k',
                                    #     lw=0.1,
                                    #     ls='--',
                                    #     headwidth=3,
                                    #     headlength=6
                                    # ),
                                    fontsize=10,
                                    ha='center',
                                    va='center'
                                    )

                    else:
                        stats_start = stats_pack
                        sm.taylor_diagram(
                            np.asarray((1, stats_pack[3])),  # mod_norm_stdev
                            np.asarray((0, stats_pack[5])),  # ncrmsd
                            np.asarray((1, stats_pack[6])),  # ccoef
                            # note change required to plot_pattern_diagram_markers.py in sm to get
                            # symbols that are letters or numbers to work
                            markersymbol=marker_dict['marker'],
                            # markercolor=modelruns_info[run]['colour'],
                            markercolor=mf_colour,
                            markersize=marker_dict['size'],
                            alpha=alpha,
                            overlay='on'
                        )
                    plt_n += 1
            SSTbuoy_scores = None

# plot the target diagrams
augment_rmsd = True
var_codes = ["T", "S"]
legend_markers_bu = []
for var_ts in var_codes:
    # determine which subplot to put plot
    if var_ts == "T":
        ax = fig.add_subplot(gs[6:11, 0:5])
        axismax = 1.5
        circles_trg = [0.5, 1.0, 1.5]
        ticks = np.asarray([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
        ax.text(-1.5, 1.7, "(c)", ha="left", va="top", fontsize=12)

    elif var_ts == 'S':
        ax = fig.add_subplot(gs[6:11, 6:11])
        axismax = 3
        circles_trg = [1.0, 2.0, 3.0]
        ticks = np.asarray([-3, -2, -1.0, 0.0, 1.0, 2.0, 3.0])
        ax.text(-2.5, 3.3, "(d)", ha="left", va="top", fontsize=12)

    overlay = 'off'
    for run in modelruns_info.keys():
        # hack
        if run == 'SalishSea1500-RUN216-altSST':
            run_sname = 'SalishSea1500-RUN216'
        else:
            run_sname = run

        if var_ts == "T":
            print("Taylor SST", run)
        elif var_ts == "S":
            print("Taylor SSS", run)
        else:
            print("No variable chosen")
            break

        path1 = modelruns_info[run]['path']
        shortcode = modelruns_info[run]['shortcode']

        # plot LH data (to do: redundant to above
        print("plotting LH data")
        LH_pic_data = 'LH_class4_' + shortcode + '_hindcast.pickle'
        with open(os.path.join(path1, LH_pic_data), 'rb') as f:
            LH_scores = pickle.load(f)

        # get stats and plot each lh
        for lh in LH_scores['obs'].keys():
            if lh not in lh_list:
                continue
            print(lh)

            # plot a marker
            try:
                marker_dict = OPTIONS_M_LH[lh]
                mf_colour = marker_dict['mf_colour']
                me_colour = marker_dict['me_colour']
            except:
                print("No marker found for " + lh + ". ")
                continue

            plt_n = 0
            for mv_avg in mv_avg_list.keys():
                print("plotting moving average " + str(mv_avg) + ' days')
                color_mvg_avg = mv_avg_list[mv_avg]['color']

                # returns buoy_bettername, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss
                stats_pack = get_target_stats_LH(lh, mv_avg, LH_scores, var_ts,
                                                 run_sname, use_nanmean, augment_rmsd)

                # if using the arrows instead of symbols then for 2nd plot use them
                if (use_arrows and plt_n > 0):
                    stats_end = stats_pack

                    x1 = stats_start[2]
                    y1 = stats_start[1]
                    x2 = stats_end[2]
                    y2 = stats_end[1]

                    ax = plt.gca()
                    if mf_colour == 'w':
                        mf_colour = '#585858'
                    ax.annotate('', xy=(x2, y2),
                                xytext=(x1, y1),
                                arrowprops=dict(
                                    arrowstyle='-|>',
                                    # head_width=0.7,
                                    # shrink_factor=0.2,
                                    mutation_scale=arrow_mutatescale,
                                    facecolor=mf_colour,
                                    edgecolor=mf_colour,
                                    lw=arrow_lw,
                                    ls='--',
                                    zorder=1
                                ),
                                # arrowprops=dict(
                                #     facecolor='k',
                                #     edgecolor='k',
                                #     lw=0.1,
                                #     ls='--',
                                #     headwidth=3,
                                #     headlength=6
                                # ),
                                fontsize=10,
                                ha='center',
                                va='center'
                                )
                else:
                    stats_start = stats_pack
                    # the sm.taylor_diagram code uses the current active axis
                    sm.target_diagram(stats_pack[1],
                                      stats_pack[2],
                                      stats_pack[3],
                                      axismax=axismax,
                                      circles=circles_trg,
                                      circlelinewidth=0.8,
                                      circlelinespec='--', # not compatible?
                                      # circlelinespec='0.75--',  # not compatible? 2024-09
                                      #circlelinespec='k--',
                                      markersymbol=marker_dict['marker'],
                                      markercolor=mf_colour,
                                      markerSize=marker_dict['size'],
                                      overlay=overlay,
                                      ticks=ticks,
                                      normalized='off',
                                      alpha=alpha)
                    if var_ts == 'T':
                        # for legend
                        mark_leg = mlines.Line2D([], [],
                                                 mfc=mf_colour,
                                                 mec='k',
                                                 marker=marker_dict['marker'],
                                                 linestyle='None',
                                                 markersize=marker_dict['size'],
                                                 label=marker_dict['label'])
                        legend_markers_bu.append(mark_leg)


                    overlay = 'on'
                plt_n += 1
        LH_scores = None

        # Plot SST from buoys
        if var_ts == 'T':
            print("plotting targets for buoys..")
            SSTbuoy_pic_data = "SST_class4_" + run + ".pickle"
            with open(os.path.join(path1, SSTbuoy_pic_data), 'rb') as f:
                SSTbuoy_scores = pickle.load(f)

            for bu in SSTbuoy_scores['obs'].keys():

                # hack
                if bu == "C46134_2001-02_2016-12.nc":
                    buoy_bettername = '46134'
                elif bu == "C46131_1992-10_2021-01.nc":
                    buoy_bettername = '46131'
                elif bu == "C46146_1992-03_2022-06.nc":
                    buoy_bettername = '46146'
                elif bu == "C46182_1989-09_1991-11.nc":
                    buoy_bettername = '46182'
                else:
                    buoy_bettername = bu

                if (buoy_bettername not in buoy_list):
                    continue
                print(buoy_bettername)

                # plot one marker
                try:
                    marker_dict = OPTIONS_M_SST[buoy_bettername]
                    mf_colour = marker_dict['mf_colour']
                    me_colour = marker_dict['me_colour']
                except:
                    print("No marker found for " + buoy_bettername + ". ")
                    continue

                plt_n = 0
                for mv_avg in mv_avg_list.keys():
                    print("plotting moving average " + str(mv_avg) + ' days')
                    color_mvg_avg = mv_avg_list[mv_avg]['color']

                    # returns buoy_bettername, stdev_obs, stdev_mod, mod_norm_stdev, crmsd, ncrmsd, ccoef, wss
                    stats_pack = get_target_stats_SST(bu, mv_avg, SSTbuoy_scores,
                                                      run_sname, use_nanmean, augment_rmsd)

                    print(stats_pack)

                    if (use_arrows and plt_n > 0):
                        stats_end = stats_pack
                        x1 = stats_start[2]
                        y1 = stats_start[1]
                        x2 = stats_end[2]
                        y2 = stats_end[1]

                        ax = plt.gca()
                        if mf_colour == 'w':
                            mf_colour = 'k'
                        ax.annotate('', xy=(x2, y2),
                                    xytext=(x1, y1),
                                    arrowprops=dict(
                                        arrowstyle='-|>',
                                        # head_width=0.7,
                                        # shrink_factor=0.2,
                                        mutation_scale=arrow_mutatescale,
                                        facecolor=mf_colour,
                                        edgecolor=me_colour,
                                        lw=arrow_lw,
                                        ls='--',
                                        zorder=1
                                    ),
                                    # arrowprops=dict(
                                    #     facecolor='k',
                                    #     edgecolor='k',
                                    #     lw=0.1,
                                    #     ls='--',
                                    #     headwidth=3,
                                    #     headlength=6
                                    # ),
                                    fontsize=10,
                                    ha='center',
                                    va='center'
                                    )

                    else:
                        stats_start = stats_pack
                        # the sm.taylor_diagram code uses the current active axis
                        sm.target_diagram(stats_pack[1],
                                          stats_pack[2],
                                          stats_pack[3],
                                          axismax=axismax,
                                          circles=circles_trg,
                                          circlelinewidth=0.8,
                                          #circlelinespec='0.75--', # not compatible?
                                          circlelinespec='k--',
                                          markersymbol=marker_dict['marker'],
                                          markercolor=mf_colour,
                                          markerSize=marker_dict['size'],
                                          overlay=overlay,
                                          ticks=ticks,
                                          normalized='off',
                                          alpha=alpha)

                        if var_ts == 'T':
                            # for legend
                            mark_leg = mlines.Line2D([], [],
                                                     mfc=mf_colour,
                                                     mec='k',
                                                     marker=marker_dict['marker'],
                                                     linestyle='None',
                                                     markersize=marker_dict['size'],
                                                     label=marker_dict['label'])
                            legend_markers_bu.append(mark_leg)
                    plt_n += 1
                    stats_pack = None
                    overlay = 'on'

# ax_tg.quiver(RMSDs_strt, Bs_strt, (RMSDs - RMSDs_strt), (Bs-Bs_strt), angles='xy', scale_units='xy', scale=1.1,
#                  color='#878787', headwidth=3, minshaft=2)

# legend
ax_leg = fig.add_subplot(gs[11:, 0:11])
ax_leg.axes.get_xaxis().set_visible(False)
ax_leg.axes.get_yaxis().set_visible(False)
ax_leg.spines['top'].set_visible(False)
ax_leg.spines['right'].set_visible(False)
ax_leg.spines['bottom'].set_visible(False)
ax_leg.spines['left'].set_visible(False)

# order the labels alphabetically then by numbers
legend_mark_sorted = sorted(legend_markers_bu, key=lambda t: t.get_label())
# for item in legend_mark_sorted:
#     print(item.get_label())
ax_leg.legend(handles = legend_mark_sorted,
                    loc="upper center",
                    ncol=4,
                    bbox_to_anchor=(0.5,1.2),
                    handletextpad=0.2,
                    fontsize=8)
plt.tight_layout()
plt.savefig('../../figs/Fig07_revised.png', dpi=300)
plt.savefig('../../figs/Fig07_revised.pdf', dpi=300)
plt.show()

