#!/bin/ksh
# in mobatexteditor ensure Unix format button is checked 
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;
#
# GO 2023-04-23
# interpolate ORAS to daily from monthly
# this script meant to be run on unix (server: graham)
# downloaded 2024-03-05
pathin="/project/6006412/goldford/data/atmos/ORAS5_JdF/"



echo interpolating to daily;

cdo -inttime,1979-01-16,12:00,1day ${pathin}votempervosaline_1m_ORAS5_1979_2018_JdF.nc ${pathin}votempervosaline_1d_ORAS5_1979_2018_JdF_interp.nc


