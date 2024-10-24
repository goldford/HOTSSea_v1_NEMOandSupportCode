#!/bin/bash
#SBATCH --job-name="plots comparison"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=1600   #4 Gb
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32 # CPUs

#module load python/3.8.10 # 2024-09 trblsht
module load StdEnv/2020
module load scipy-stack/2020a
#source pyenv/bin/activate # 2024-09 trblsht
#python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml
#python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG -j 32
# CTD SST MCTD LH TG
#module load scipy-stack/2023b
#python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG
#python gitlab_pyap_fork20230331/bin/scores.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG -j 3

#python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/compare_several_GO.yaml -i CTD -j 3
python gitlab_pyap_fork20230331/bin/plotcomparisons_GO.py gitlab_pyap_fork20230331/config/mid002/SS1500/compare_several_GO.yaml -i CTD -j 3 # not sure why above not working, trying this 2024-09

