import datetime
import numpy as np
import os
import scipy.signal as signal

from analysispkg import pkg_data
from analysispkg import pkg_obs
from analysispkg import pkg_signal
from analysispkg import pkg_statistics
from analysispkg import pkg_utils


def analyze_SST(opt):

    casename = opt['casename']
    print('Processing data for %s...' % casename)

    outpath = os.path.join(opt['dir_process'], 'SST')
    os.makedirs(outpath, exist_ok=True)

    for period, (ana_start, ana_end) in opt['analysis']['periods'].items():
        print("Working on period {}".format(period))

        print('Loading obs database...')
        inst = make_instrument_dict(opt, ana_start, ana_end)

        print('Reading model data...')
        dat = {}
        dat[casename] = read_model(opt, inst, ana_start, ana_end)

        print('Amalgamating data from all sources, and filtering/interpolating for fair comparison...')
        class4 = make_class4().run(dat, inst, opt)

        #filename = 'SST_class4_' + opt['runid'] + "_" + period + ".pickle"
        outpath, filename = pkg_utils.get_process_path(opt, 'SST', period)
        os.makedirs (outpath, exist_ok=True)

        print('Saving {} ...'.format(filename))
        #pkg_data.save_pickle_and_mat(os.path.join(outpath, filename), class4)
        pkg_data.save_pickle_and_mat(opt, filename, class4)

    print('Done!')


def make_instrument_dict(opt, ana_start, ana_end):
    ostart = datetime.datetime(2030, 12, 31)
    oend = datetime.datetime(1900, 1, 1)

    instr_type = 'SST'

    # Load observations index and filter by location, time, and exclude list
    observations = pkg_obs.observations_load_and_filter(opt, instr_type, ana_start, ana_end)

    inst = {}

    for ob in observations:
        obsfile = ob['filename']
        print (ob)
        obs = pkg_obs.load_observation(obsfile, 'SST')
        obs = pkg_obs.truncate_obs_timeseries(obs, ana_start, ana_end)
        obs['time'], obs['temperature'] = pkg_signal.insert_nans_in_gaps(obs['time'], obs['temperature'])

        # NaN removal?
        #     c2 = np.where(np.isnan(obs['T']))[0]
        #     obs['time'] = np.delete(obs['time'], c2)
        #     obs['T'] = np.delete(obs['T'], c2)

        if obs['time'][0] < ostart:
            ostart = obs['time'][0]
        if obs['time'][-1] > oend:
            oend = obs['time'][-1]

        _, fname = os.path.split(obsfile)
        station = fname.split('_SST')[0]

        write = {'lat': obs['lat'], 'lon': obs['lon'], 'start': obs['time'][0], 'end': obs['time'][-1], 'obs': obs,
                 'obs_file(s)': [obsfile], 'shortname':obs.get('shortname', None)}

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

    obs = inst[s]['obs']
    moddata, pfile = pkg_data.load_extraction(opt, obs)
    moddata = pkg_obs.truncate_mod_timeseries(moddata, ana_start, ana_end, ['time', 'temperature'])
    if moddata is not None:
        res[s] = moddata
    return res


class make_class4():
    def filt_var(self, var, src_dt, cut_time, butter_order=2):
        fs = 1. / src_dt
        fc = 1. / cut_time
        w = fc / (fs / 2.)
        if w >= 1.:
            return var
        else:
            b, a = signal.butter(butter_order, w)
            filt = signal.filtfilt(b, a, var)

        return filt

    def compare_timesteps(self, dat, s):
        dtmax = 0.
        max_src = None
        dts = {}
        tend = datetime.datetime(2030, 12, 31)
        tstart = datetime.datetime(1900, 1, 1)
        iterator = list(dat.keys())
        for src in iterator:
            if s in dat[src].keys():
                if dat[src][s]['time'][-1] < tend:
                    tend = dat[src][s]['time'][-1]
                if dat[src][s]['time'][0] > tstart:
                    tstart = dat[src][s]['time'][0]
                dta = np.array([d.total_seconds() for d in np.diff(dat[src][s]['time'])])
                if src == 'obs' and np.std(dta) > 60.:
                    print(
                        'WARNING: observations are at inconsistent time intervals for %s. Linearly interpolating to even time interval.' % s)
                    tRaw = np.array([(t - dat[src][s]['time'][0]).total_seconds() for t in dat[src][s]['time']])
                    tI = np.arange(0., (dat[src][s]['time'][-1] - dat[src][s]['time'][0]).total_seconds(),
                                   np.mean(dta[dta != dta.max()]))
                    dat[src][s]['temperature'] = np.interp(tI, tRaw, dat[src][s]['temperature'])
                    dat[src][s]['time'] = np.array([dat[src][s]['time'][0] + datetime.timedelta(seconds=t) for t in tI])
                    dts[src] = np.mean(dta[dta != dta.max()])
                else:
                    dts[src] = np.mean(dta)
                if dts[src] > dtmax:
                    dtmax = dts[src]
                    max_src = src

        return dts, dtmax, max_src, tstart, tend

    def filter_and_interpolate(self, dat, int_filt_dat, s):
        dts, dtmax, max_src, tstart, tend = make_class4().compare_timesteps(dat, s)

        for src in dat.keys():
            if s in dat[src].keys():
                if src not in int_filt_dat.keys():
                    int_filt_dat[src] = {}
                if s not in int_filt_dat[src].keys():
                    int_filt_dat[src][s] = {}

        tI = np.array([(t - tstart).total_seconds() for t in dat[max_src][s]['time']])
        if s in dat[max_src].keys():
            int_filt_dat[max_src][s]['temperature'] = dat[max_src][s]['temperature']
        for src in dat.keys():
            if s in dat[src].keys():
                if 'time' not in int_filt_dat[src][s].keys():
                    int_filt_dat[src][s]['time'] = dat[max_src][s]['time']
                if src != max_src:
                    tS = np.array([(t - tstart).total_seconds() for t in dat[src][s]['time']])
                    nanflag = False
                    var = dat[src][s]['temperature']
                    c = np.where(np.logical_or(np.isnan(var), var == 0.))[0]
                    if len(c) == len(var):
                        print('All NaNs for %s' % s)
                        nanflag = True
                        if nanflag:
                            print('Couldn\'t find any useful data from %s for %s. Deleting this entry from results.' % (
                                src, s))
                            del int_filt_dat[src][s]
                    else:
                        if len(c) > 0:
                            tS = np.delete(tS, c)
                            var = np.delete(var, c)
                        if dts[src] < dtmax:
                            var = make_class4().filt_var(var, dts[src], dtmax)
                        tI_nan_mask = np.zeros_like(tI).astype(int)
                        for i,t in enumerate(tI):
                            if np.abs(tS - t).min() > 86400.:
                                tI_nan_mask[i] = 1
                        int_filt_dat[src][s]['temperature'] = np.interp(tI, tS, var)
                        int_filt_dat[src][s]['temperature'][tI_nan_mask == 1] = np.nan

        # calculate model scores
        for src in int_filt_dat.keys():
            if src != 'obs':
                for s in int_filt_dat[src].keys():
                    if s in int_filt_dat['obs'].keys():  # some stations can be missing from obs series, the ones with no valid data
                        scores = pkg_statistics.residual_stats(int_filt_dat['obs'][s]['temperature'],
                                                               int_filt_dat[src][s]['temperature'])
                        int_filt_dat[src][s]['scores'] = scores

        return int_filt_dat

    def run(self, dat, inst, opt):
        dat['obs'] = {}
        int_filt_dat = {}
        for s in inst.keys():
            if s not in ['start', 'end']:
                dat['obs'][s] = inst[s]['obs']
                int_filt_dat = make_class4().filter_and_interpolate(dat, int_filt_dat, s)

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
