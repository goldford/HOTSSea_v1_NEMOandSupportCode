import argparse

from analysispkg import pkg_scan
from analysispkg import pkg_utils


def main():

    parser = argparse.ArgumentParser()
    # Required positional parameters
    parser.add_argument("optfile", type=str,
                        help="YAML file containing parameters.")
    parser.add_argument("-j", type=pkg_utils.positive_int,
                        help="Number of processors to use for parallelized calculation.")
    parser.add_argument("-f", action='store_true',
                        help="Force rescan (discard any existing index), default is false.")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity level (-v, -vv or -vvv for 1, 2 or 3)")
    args = parser.parse_args()

    # Load parameters from YAML file
    opt = pkg_utils.load_options(args.optfile)
    if args.verbose > 0:
        opt['verbose'] = args.verbose

    if args.j is not None:
        opt['nproc'] = args.j
        opt['parallel'] = opt['nproc'] > 1

    pkg_utils.prologue(opt)

    pkg_scan.scan(opt,force=args.f)

    print("scan.py complete")

