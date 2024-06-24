import argparse

from analysispkg import pkg_csv
from analysispkg import pkg_utils


def main():
    parser = argparse.ArgumentParser()
    # Required positional parameters
    parser.add_argument("optfile", type=str,
                        help="YAML file containing parameters.")
    # Optional parameters
    parser.add_argument("--instruments", "-i", type=str, nargs='+', choices=pkg_utils.INSTRUMENTS,
                        help="Scores for specific instruments. Default is to score instruments specifies in config YAMLs.")
    parser.add_argument("-j", type=pkg_utils.positive_int, default=1,
                        help="Number of processors to use for parallelized calculation. Currently the scoring is not parallelized, so this should be set to 1.")
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

    if args.instruments is not None:
        opt['instruments'] = args.instruments

    pkg_utils.prologue(opt)

    for instrument in opt['instruments']:
        print("Scores for instrument {}".format(instrument))

        if instrument in ['CTD']:
            pkg_csv.write_scores_ctd(opt)

        if instrument in ['LH']:
            pkg_csv.write_scores_lh(opt)

        if instrument in ['SST']:
            pkg_csv.write_scores_sst(opt)

    print("scores.py complete.")
