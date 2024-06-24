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
                        help="Analyze specific instruments. Default is to analyze instruments specified via config YAML.")
    parser.add_argument("-j", type=pkg_utils.positive_int,
                        help="Number of processors to use for parallelized calculation.")
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help="Launch pdb on unhandled exception")
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

    if args.debug:
        pkg_utils.debug_switch[0]=True

    # Passed args validation. If we provided "--submit", then we submit it to a node
    if len(sysnames) > 0:
        if args.submit is True:
            myargs = sys.argv.copy()
            myargs.remove('--submit')
            logdir = opt['dir_logs']
            sysname = args.c
            pkg_submit.submit(sysname,myargs,logdir,nodes=1)
            sys.exit()

    if args.instruments is not None:
        opt['instruments'] = args.instruments

    pkg_utils.prologue(opt)

    analyze(opt)

    print("analyze.py complete")


@pkg_utils.debug
def analyze(opt):
    for instrument in opt['instruments']:
        print("Analyze instrument {}".format(instrument))

        if instrument in ['CTD']:
            from analysispkg import pkg_analyze_cast
            pkg_analyze_cast.analyze_casts(opt)

        if instrument in ['LH']:
            from analysispkg import pkg_analyze_LH
            pkg_analyze_LH.analyze_LH(opt)
        if instrument in ['SST']:
            from analysispkg import pkg_analyze_SST
            pkg_analyze_SST.analyze_SST(opt)
