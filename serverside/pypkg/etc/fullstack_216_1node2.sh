#!/bin/bash

#module load StdEnv/2020
#module load scipy-stack/2020a
#python gitlab_pyap_fork20230331/bin/scan.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml
python gitlab_pyap_fork20230331/bin/extract.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG -j 8
# CTD SST MCTD LH TG
module load scipy-stack/2020b
python gitlab_pyap_fork20230331/bin/analyze.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG
#python gitlab_pyap_fork20230331/bin/scores.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG -j 3
#python gitlab_pyap_fork20230331/bin/plots.py gitlab_pyap_fork20230331/config/mid002/SS1500/SalishSea1500-RUN216.yaml -i CTD SST MCTD LH TG -j 3

