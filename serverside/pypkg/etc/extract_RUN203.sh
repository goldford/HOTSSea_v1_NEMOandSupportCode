#!/bin/bash
#SBATCH --job-name="ctd extract pyap test"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=2:00:00
#SBATCH --ntasks=3 # CPUs -GO

python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203.yaml -i TG -j 3
