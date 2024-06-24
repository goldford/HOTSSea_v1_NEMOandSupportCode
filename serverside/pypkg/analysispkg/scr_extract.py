import argparse
import sys

from analysispkg import pkg_submit
from analysispkg import pkg_utils


def main():
    parser = argparse.ArgumentParser()
    # Required positional parameters
    parser.add_argument("optfile", type=str,
                        help="YAML file containing parameters.")
    # Optional parameters
    parser.add_argument("--instruments", "-i", type=str, nargs='+', choices=pkg_utils.INSTRUMENTS,
                        help="Extract specified instrument types. Default is to extract instruments specified via YAMLs.")
    parser.add_argument("-j", type=pkg_utils.positive_int,
                        help="Number of processors to use for parallelized extraction.")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity level (-v, -vv or -vvv for 1, 2 or 3)")
    parser.add_argument("-t", "--test", action='store_true',
                        help="Test mode (run with truncated observations for rapid testing).")
    sysnames = pkg_submit.get_sysnames()
    if len(sysnames) > 0:
        parser.add_argument("--submit", action='store_true',
                            help="Submit to the queue using one full node.")
        parser.add_argument('-c', choices=sysnames, default=sysnames[0],
                            help="Select which system to submit to; default {}.".format(sysnames[0]))
    args = parser.parse_args()

    # Load parameters from YAML file
    opt = pkg_utils.load_options(args.optfile)
    if args.verbose > 0:
        opt['verbose'] = args.verbose
    opt['test'] = args.test

    if args.j is not None:
        opt['nproc'] = args.j
        opt['parallel'] = opt['nproc'] > 1

    # Passed args validation. If we provided "--submit", then we submit it to a node
    if len(sysnames) > 0:
        if args.submit is True:
            sysname = args.c
            myargs = sys.argv.copy()
            myargs.remove('--submit')
            logdir = opt['dir_logs']
            pkg_submit.submit(sysname,myargs,logdir,nodes=1)
            sys.exit()

    if args.instruments is not None:
        opt['instruments'] = args.instruments

    pkg_utils.prologue(opt)

    for instrument in opt['instruments']:
        print("Extracting instrument {}".format(instrument))

        if instrument in ['CTD']:
            from analysispkg import pkg_extract_cast
            pkg_extract_cast.extract_model(opt)

        if instrument in ['SST','LH']:
            from analysispkg import pkg_extract_timeseries
            pkg_extract_timeseries.extract_model_timeseries(opt,instrument)
    print("extract.py complete")
