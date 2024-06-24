#!/bin/bash
#SBATCH --job-name="all process stack"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=04:00:00
#SBATCH --ntasks=4 # CPUs

module load scipy-stack/2020a
python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_xcheck.yaml
python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_xcheck.yaml -i CTD -j 4
# CTD SST MCTD LH TG
module load scipy-stack/2020b
python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_xcheck.yaml -i CTD -j 4
python gitlab_pyap_fork20230331/bin/scores.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_xcheck.yaml -i CTD -j 4
python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN203_xcheck.yaml -i CTD -j 4

