#!/bin/bash
#SBATCH --job-name="ctd analyze"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=01:00:00
#SBATCH --ntasks=10 # CPUs


python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_GO.yaml -i CTD -j 10
python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD -j 10
