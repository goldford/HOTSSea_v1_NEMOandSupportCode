import csv
import os
import pickle

import numpy as np
from matplotlib.path import Path as Path

from analysispkg import pkg_analyze_cast, pkg_utils

csvopts = {'delimiter': ',',
           'quotechar': '|',
           'lineterminator': '\n'}


def write_scores_sst(opt):
    """ Write a csv file with statistical scores for SST data.
    """
    typ = 'SST'
    for period in opt['analysis']['periods'].keys():
        _, filename = pkg_utils.get_process_path(opt, 'SST', period)
        with open(filename, 'rb') as fid:
            fulldata = pickle.load(fid)
        if len(fulldata) == 0:
            continue

        # the only model among keys; the other key is 'obs'
        model_id = (set(fulldata.keys()) - {'obs'}).pop()
        fulldata = fulldata[model_id]

        # Ensure output has a place to go
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        os.makedirs(out_dir, exist_ok=True)

        for statistic in opt['score']['statistics']:
            codes, locations, stats = [], [], []
            for instrument, data in fulldata.items():  # stations
                if 'scores' in data['filt_interp_data'].keys():
                    codes.append(instrument)
                    locations.append([data['lon'], data['lat']])
                    stats.append(data['filt_interp_data']['scores'][statistic])

            # Output in csv file
            fn = os.path.join(out_dir, pkg_utils.csv_file(typ, period, 'scores_temperature', statistic, opt))
            with open(fn, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile, **csvopts)
                csvwriter.writerow(['Station name', 'Obs. Lon', 'Obs. Lat', statistic])
                for stn, loc, stat in zip(codes, locations, stats):
                    csvwriter.writerow([stn, loc[0], loc[1], np.round(stat, 3)])

    print('pkg_csv.write_scores_sst complete.')


def write_scores_lh(opt):
    """ Write a csv file with statistical scores for LH data.
    TODO: see if we can converge this code with SST
    """
    typ = 'LH'
    for period in opt['analysis']['periods'].keys():
        file_summary = os.path.join(opt['dir_process'], typ,
                                    typ + '_class4_' + opt['runid'] + '_' + period + '.pickle')
        with open(file_summary, 'rb') as fid:
            fulldata = pickle.load(fid)

        # If no data, no csv should be written.
        if fulldata:
            obsdata = fulldata['obs']
            # the only model among keys; the other key is 'obs'
            try:
                model_id = (set(fulldata.keys()) - {'obs'}).pop()
            except:
                print ('    No valid model data for period {}'.format(period))
                continue
            else:
                fulldata = fulldata[model_id]

                # Ensure output has a place to go
                out_dir = pkg_utils.get_plot_dir(opt, typ, period)
                os.makedirs(out_dir, exist_ok=True)

                for var in ['temperature', 'salinity']:
                    for statistic in opt['score']['statistics']:
                        codes, locations, stats = [], [], []
                        for instrument in fulldata.keys():  # stations
                            if instrument in obsdata.keys():  #TODO These series should be skipped from model dict. Then this check becomes obsolete.
                                data = fulldata[instrument]
                                codes.append(instrument)
                                locations.append([data['lon'], data['lat']])
                                stats.append(data['filt_interp_data']['scores'][var][statistic])

                        # Output in csv file
                        fn = os.path.join(out_dir, pkg_utils.csv_file(typ, period, 'scores_' + var, statistic, opt))
                        with open(fn, 'w', newline='') as csvfile:
                            csvwriter = csv.writer(csvfile, **csvopts)
                            csvwriter.writerow(['Station name', 'Obs. Lon', 'Obs. Lat', statistic])
                            for stn, loc, stat in zip(codes, locations, stats):
                                csvwriter.writerow([stn, loc[0], loc[1], np.round(stat, 3)])
        else:
            print ('    No LHs for analysis period {}'.format(period))

    print('pkg_csv.write_scores_lh complete.')

def write_scores_ctd(opt):
    """ Output statistical scores to a csv file.

    Parameters
    ----------
    opt
    """
    domain_file = opt['analysis']['CTD']['domain_file']
    _, regcoords = pkg_analyze_cast.read_CTD_domain(domain_file)

    typ = 'CTD'
    for period in opt['analysis']['periods'].keys():
        print("Working on scores for {}".format(period))
        file_summary = os.path.join(opt['dir_process'], typ,
                                    'CTDcast_metrics_' + period + '.pickle')
        with open(file_summary, 'rb') as fid:
            pdata = pickle.load(fid)
        if pdata is None:
            continue
        data = pdata[opt['casename']]

        # Ensure output has a place to go
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        os.makedirs(out_dir, exist_ok=True)

        scores = data['scores']
        if len(scores) == 0:
            print("No CTD scores for period {}, moving on".format(period))
            continue

        metrics = set()
        for m in scores:
            if scores[m] is not None:
                metrics = metrics.union(scores[m].keys())  # available metrics
        # select available metric from those that are requested
        requested_metrics = [m for m in opt['score']['statistics'] if m in metrics]
        for var in ['T', 'S']:
            for statistic in requested_metrics:  # for each requested available metric
                # Collect stats for all stations
                codes, locations, stats = [], [], []
                for stn in scores.keys():
                    #TODO (2020-10-01): Why Nones appear in scores? e.g. for NQ1_CastCTD_2017-03-19_1316h
                    if scores[stn] is not None:
                        codes.append(stn)

                        locations.append(scores[stn]['ModelLocation'])

                        # get stats for the entire profile, 'Full';
                        # also available in standard depth ranges, 0-10, 10-30, ...
                        stats.append(scores[stn][statistic]['Full'][var])

                pts = np.array(locations)

                # Output in csv file
                filename = pkg_utils.csv_file(typ, period, 'scores', statistic + '_' + var, opt)
                filepath = os.path.join(out_dir, filename)
                with open(filepath, 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, **csvopts)
                    csvwriter.writerow(['Station name', 'Obs. Lon', 'Obs. Lat', 'Region', statistic])
                    for region in regcoords.keys():
                        flags = Path(regcoords[region]).contains_points(pts) # find casts in this region
                        for i, (stn, loc, stat) in enumerate(zip(codes, locations, stats)):
                            if flags[i]:
                                csvwriter.writerow([stn, loc[1], loc[0], region, np.round(stat, 3)])
        print("pkg_csv.write_scores_ctd complete for period {}".format(period))

