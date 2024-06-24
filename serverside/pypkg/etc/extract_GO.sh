#!/bin/bash
#SBATCH --job-name="ctd rxtrct"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=02:00:00
#SBATCH --ntasks=3 # guess CPUs? -GO

python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j 3
