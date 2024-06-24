#!/bin/bash
#SBATCH --job-name="all process stack"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=1200   #4 Gb
#SBATCH --time=07:40:00 # needs 4-6 hrs for extract j 10
#SBATCH --ntasks=10 # CPUs

module load StdEnv/2020
module load scipy-stack/2020a
#python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml
python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j  10 -v
# CTD SST MCTD LH TG CM

module load scipy-stack/2020b
#python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j 10
#python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j 10
#python gitlab_pyap_fork20230331/bin/scores.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j 10
#python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN100_GO.yaml -i CM -j 10