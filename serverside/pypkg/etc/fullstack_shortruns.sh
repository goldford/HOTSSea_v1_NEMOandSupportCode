#!/bin/bash
#SBATCH --job-name="ctd stack"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=02:00:00
#SBATCH --ntasks=4 # CPUs

#python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml 
python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CTD SST LH TG -j 4
python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CTD SST LH TG -j 4
python gitlab_pyap_fork20230331/bin/scores.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CTD SST LH TG -j 4
python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CTD SST LH TG -j 4

