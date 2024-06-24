import os
import sys
import glob
import shutil
import datetime
import numpy as np
import netCDF4 as nc
from analysispkg import pkg_utils, pkg_obs, pkg_data


ZERO_DELTA = datetime.timedelta(0)
GARBAGE_SUBDIR = 'garbage_update_archive'
LOGS_SUBDIR = 'logs_update_archive'


def update_instrument(archive, instrument, updated_files_dir,
                      dry_run=False, out_subdir='', dist_threshold_km=0.1, inst_dep_threshold=5):
    """ Update data files in the obs archive from a directory with new files.

    Copies new files to the archive.
    Compares new and archived files by location, time, and instrument depth (if applicable)
    and moves old archived files that are overlapped by the new ones to a garbage dir,
    ARCHIVE_ROOT/garbage_update_archive.
    An archive file that has ANY overlap in time with a new file is considered a match
    and is moved to the garbage.
    Changes are saved in a log file in ARCHIVE_ROOT/logs_update_archive.

    The code should be run from the account owning the data archive, i.e. sdfo600.
    So far has been run with 'CTD'.
    TODO: Test with other instruments.

    Uses brute force search to match data.

    Parameters
    ----------
    archive : {'w' | 'e' | 'test'}
        Data archive.
    instrument : str, {'CTD' | ...}
        Instrument type.
    updated_files_dir : str
        Dir with new files.
    dry_run : bool, optional
        Does not do any file operations, only shows what will be done. Default is `False`.
    out_subdir : str, optional
        Subdirectory for the new files in the archive,
        relative to archive/instrument dir. Default is no subdirectory.
    dist_threshold_km : float, optional
        Threshold for location matching, km. Default is 0.1 km.
    inst_dep_threshold : float, optional
        Threshold for instrument depth matching, m. Default is 5 m.
    """
    # Get archive config with root dir for the archive
    configpath = os.path.normpath(os.path.dirname(__file__) + "/../config")
    if archive.lower() in ['w', 'west']:
        opt = pkg_utils.load_yaml(os.path.join(configpath, 'obs_west.yaml'))
    else:
        raise RuntimeError('Unknown archive.')

    # Get archive dir for the instrument, set up garbage and log dirs
    instr_dir = os.path.join(opt['root'], instrument)
    out_dir = os.path.join(instr_dir, out_subdir)
    garbage_dir = os.path.join(opt['root'], GARBAGE_SUBDIR, instrument)
    logs_dir = os.path.join(opt['root'], LOGS_SUBDIR)
    if not dry_run:
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(garbage_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)

    # Load index with existing obs data
    ind = pkg_utils.load_index(opt, instrument)
    alon = [k['lon'] for k in ind]
    alat = [k['lat'] for k in ind]
    ind_arr = np.array(ind)

    # Go over new files
    # NOTE: Searches all subdirs in `updated_files_dir`,
    # but then copies to archive without subdirs,
    # so file names in `updated_files_dir` should be unique.
    search_string = os.path.join(updated_files_dir, "**/", "*.nc")
    files = glob.glob(search_string, recursive=True)
    new_files = []  # to hold data files with no matches in the archive
    arch_to_garbage = []  # to hold archived files that are superseded by new ones
    arch_supersede = []  # to hold new files that supersede each file in arch_to_garbage
    supersede_summary = []  # to hold summaries for arch_supersede
    skipped = []  # to hold skipped new files, e.g. no data
    for f in files:
        # print(f)
        station, lon, lat, inst_depth, tstart, tend = file_meta(f, instrument)
        if tstart is None:
            skipped.append(f)
            continue
        coord_match, coord_delta = match_location(lon, lat, alon, alat, dist_threshold_km)
        ind_coord_matched = ind_arr[coord_match]
        coord_delta = coord_delta[coord_match]
        # if len(ind_coord_matched) > 0:
        #     print(f)
        #     print(len(ind_coord_matched))
        found_match = False
        # go over archived files
        for indf, coord_delta1 in zip(ind_coord_matched, coord_delta):
            # print('match_time')
            overlap, delta_arc = match_time(tstart, tend, indf)
            # print('   ', overlap, delta_arc)
            if overlap > ZERO_DELTA:
                # match instrument depth
                dep_match, dep_delta = True, 0
                if dep_match:
                    # loc, time, instr dep all match
                    found_match = True
                    overlap_perc = round(100 * overlap.total_seconds() / delta_arc.total_seconds(), 2)
                    # save summary
                    summary = 'overlaps archived by: {} ({}%)'.format(overlap, overlap_perc)
                    if coord_delta1 > 0:
                        summary += '; coord delta: {}'.format(coord_delta1)
                    if dep_delta > 0:
                        summary += '; instr depth delta: {}'.format(dep_delta)
                    arch_to_garbage.append(indf['filename'])
                    arch_supersede.append(f)
                    supersede_summary.append(summary)

        # if not found_match:
        new_files.append(f)

    # summary of file changes
    if dry_run:
        sum_out = sys.stderr
    else:
        log_name = datetime.datetime.utcnow().strftime("%Y-%m-%d_%H-%M-%S") + '.log'
        sum_out = open(os.path.join(logs_dir, log_name), 'w')

    def print_ln(*args):
        print(*args, file=sum_out)

    print_ln('==============================')
    print_ln('REPLACED FILES in the archive:')
    print_ln('==============================')
    if len(arch_to_garbage) == 0:
        print_ln('None')
    else:
        arch_to_garbage_unique, ui = np.unique(arch_to_garbage, return_inverse=True)
        arch_supersede = np.array(arch_supersede)
        supersede_summary = np.array(supersede_summary)
        for k, f1 in enumerate(arch_to_garbage_unique):
            if not dry_run:
                f1save = f1.replace(instr_dir, garbage_dir)
                os.makedirs(os.path.split(f1save)[0], exist_ok=True)  # ensure such subdir exists in garbage
                shutil.move(f1, f1save)
            print_ln(f1)
            print_ln('    Superseded by:')
            sup = arch_supersede[ui == k]
            supsum = supersede_summary[ui == k]
            isort = np.argsort(sup)
            for sup1, sum1 in zip(sup[isort], supsum[isort]):
                print_ln('   ', sup1)
                print_ln('      ', sum1)

    print_ln('\n')
    print_ln('================================')
    print_ln('NEW FILES copied to the archive:')
    print_ln('================================')
    if len(new_files) == 0:
        print_ln('None')
    else:
        new_files.sort()
        for f1 in new_files:
            if not dry_run:
                shutil.copy2(f1, out_dir)
            print_ln(f1)


    if len(skipped) > 0:
        print_ln('\n')
        print_ln('============================')
        print_ln('SKIPPED NEW FILES (no data):')
        print_ln('============================')
        skipped.sort()
        for f1 in skipped:
            print_ln(f1)
    if not dry_run:
        sum_out.close()
    # return new_files, arch_to_garbage, arch_supersede


def match_location(lon, lat, alon, alat, dist_threshold_km):
    # match location lon, lat with a set of points alon, alat
    coord_delta = pkg_data.haversine(lon, lat, alon, alat)
    coord_delta[coord_delta < 0.001] = 0  # ignore diffs < 1 m
    coord_match = coord_delta < dist_threshold_km
    return coord_match, coord_delta


def match_time(tstart, tend, indf):
    # match time
    delta_arc = indf['time'][-1] - indf['time'][0]
    delta_st = tstart - indf['time'][0]
    delta_en = indf['time'][-1] - tend
    if delta_st > delta_arc or delta_en > delta_arc:
        overlap = ZERO_DELTA
    else:
        if delta_st < ZERO_DELTA:
            delta_st = ZERO_DELTA
        if delta_en < ZERO_DELTA:
            delta_en = ZERO_DELTA
        overlap = delta_arc - delta_st - delta_en
    return overlap, delta_arc


def match_depth(inst_depth, indf, inst_dep_threshold):
    with nc.Dataset(indf['filename']) as ncf:
        idep_arch = ncf.variables['inst_depth'][:]
    dep_delta = abs(idep_arch - inst_depth)
    dep_match = dep_delta < inst_dep_threshold
    return dep_match, dep_delta


def file_meta(f, instrument):
    station = os.path.basename(f).split('_')[0]
    with nc.Dataset(f) as ncf:
        lon = ncf.variables["longitude"][:]
        lat = ncf.variables["latitude"][:]
        # print(ncf)
        inst_depth = None
        tstart, tend = read_time(ncf)
    return station, lon, lat, inst_depth, tstart, tend


def read_time(ncid, timename='time', calendar='standard'):
    """Reads the time vector from a netcdf file and converts it to a datetime object.
    Arguments:
    ncid: An opened netcdf file handle
    timename: name of the time vector; default 'time'
    calendar: If calendar is not included in the netcdf attributes, use this value (default standard)
    """
    if timename not in ncid.variables:
        raise RuntimeError('Time variable {} not found.'.format(timename))

    units = ncid.variables[timename].units
    try:
        cal = ncid.variables[timename].calendar
    except:
        cal = calendar

    # convert the numbers from the file to datetime objects
    time_var = ncid.variables[timename]
    if len(time_var) == 0:
        tstart, tend = None, None
    else:
        t = nc.num2date(time_var[0], units, cal)
        tstart = datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second)
        t = nc.num2date(time_var[-1], units, cal)
        tend = datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second)
    return tstart, tend
