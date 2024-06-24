import datetime
import os

import numpy as np

from analysispkg import pkg_data, pkg_obs, pkg_statistics


# TODO: try and consolidate / converge some of this code with SST analysis

def analyze_LH(opt):

    casename = opt['casename']
    print('Processing data for %s...' % casename)

    outpath = os.path.join(opt['dir_process'], 'LH')
    os.makedirs(outpath, exist_ok=True)

    for period, (ana_start, ana_end) in opt['analysis']['periods'].items():
        print("Working on period {}".format(period))

        print('Loading obs database...')
        inst = make_instrument_dict(opt, ana_start, ana_end)

        print('Reading model data...')
        dat = {}
        dat[casename] = read_model(opt, inst, ana_start, ana_end)

        print('Amalgamating data from all sources, and filtering/interpolating for fair comparison...')
        class4 = make_class4().run(dat, inst)

        filename = 'LH_class4_' + opt['runid'] + "_" + period + ".pickle"
        print('Saving {} ...'.format(filename))
        pkg_data.save_pickle_and_mat(opt, os.path.join(outpath, filename), class4)

    print('Done!')


def make_instrument_dict(opt, ana_start, ana_end):
    ostart = datetime.datetime(2030, 12, 31)
    oend = datetime.datetime(1900, 1, 1)

    # Load observations index and filter by location, time, and exclude list
    observations = pkg_obs.observations_load_and_filter(opt, 'LH', ana_start, ana_end)

    inst = {}

    for ob in observations:
        obsfile = ob['filename']
        obs = pkg_obs.load_observation(obsfile, 'LH')
        if obs is None:
            continue
        obs = pkg_obs.truncate_obs_timeseries(obs, ana_start, ana_end)
        if obs is None:
            continue

        if obs['time'][0] < ostart:
            ostart = obs['time'][0]
        if obs['time'][-1] > oend:
            oend = obs['time'][-1]

        _, fname = os.path.split(obsfile)
        station = fname.split('_SST')[0]

        write = {'lat': obs['lat'], 'lon': obs['lon'], 'start': obs['time'][0], 'end': obs['time'][-1], 'obs': obs,
                 'obs_file(s)': [obsfile]}

        if station not in inst.keys():
            inst[station] = write

    inst['start'] = ostart
    inst['end'] = oend

    return inst


def arrange_results(results):
    full_results = {}
    for res in results:
        for s in res.keys():
            full_results[s] = res[s]

    return full_results


def read_model(opt,inst, ana_start, ana_end):
    # read series for all stations and depths
    results = []
    for s in inst.keys():  # for each station
        if s not in ['start','end']:
            for z in inst[s].keys():
                results.append( read_one_instrument(opt,inst,s, ana_start, ana_end) )
    # arrange in list of lists by station
    out = arrange_results(results)
    return out


def read_one_instrument(opt, inst, s, ana_start, ana_end):
    res = {}
    res[s] = {}

    obs = inst[s]['obs']
    moddata, pfile = pkg_data.load_extraction(opt, obs)
    moddata = pkg_obs.truncate_mod_timeseries(moddata, ana_start, ana_end, ['time', 'temperature', 'salinity'])
    if moddata is not None:
        res[s] = moddata
    return res


class make_class4():
    def run(self, dat, inst):
        dat['obs'] = {}

        # Build dict to hold results
        int_filt_dat = {}
        for src in dat.keys():
            int_filt_dat[src] = {}
            for s in dat[src].keys():
                int_filt_dat[src][s] = {}
        int_filt_dat['obs'] = {}

        # Loop over stations
        stations = set(inst.keys()) - set(['start', 'end'])
        for s in stations:
            dat['obs'][s] = inst[s]['obs']
            int_filt_dat['obs'][s] = inst[s]['obs']
            # calculate model scores
            models = set(dat.keys()) - set(['obs'])
            for src in models:
                int_filt_dat[src][s]['scores'] = {}
                for tsname in ['temperature','salinity']:
                    tobs = dat['obs'][s]['time']
                    print("GO check.")
                    print(dat[src][s].keys())
                    tmod = dat[src][s]['time']
                    t0 = min(np.min(tobs), np.min(tmod))
                    ttobs = np.array([(t - t0).total_seconds() for t in tobs])
                    ttmod = np.array([(t - t0).total_seconds() for t in tmod])
                    int_filt_dat[src][s][tsname]=np.interp(ttobs, ttmod, dat[src][s][tsname])
                    int_filt_dat[src][s]['time'] = tobs
                    scores = pkg_statistics.residual_stats(dat['obs'][s][tsname],
                                                           int_filt_dat[src][s][tsname])
                    int_filt_dat[src][s]['scores'][tsname] = scores

    # Put everything in class4 dict
        class4 = {}
        for src in int_filt_dat.keys():  # models
            class4[src] = {}
            for s in int_filt_dat[src].keys():  # stations
                class4[src][s] = {}
                for v in inst[s].keys():  # variables
                    if v == 'obs':
                        continue
                    else:
                        class4[src][s][v] = inst[s][v]
                class4[src][s]['raw_data'] = dat[src][s]
                class4[src][s]['filt_interp_data'] = int_filt_dat[src][s]

        return class4
