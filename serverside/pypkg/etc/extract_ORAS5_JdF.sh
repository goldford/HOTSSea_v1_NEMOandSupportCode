#!/bin/bash
#SBATCH --job-name="ctd extract pyap ORAS5"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=12 # CPUs

python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/ORAS5_JDF_Eval.yaml -i CTD -j 12
