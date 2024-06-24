#!/bin/bash
#SBATCH --job-name="ctd extract pyap test"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=02:00:00
#SBATCH --ntasks=3 # CPUs

python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/ORAS5_SalishSea_Eval.yaml -j 10
