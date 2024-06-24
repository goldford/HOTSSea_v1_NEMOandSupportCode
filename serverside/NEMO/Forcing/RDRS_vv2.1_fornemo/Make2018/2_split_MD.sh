#!/bin/bash
module load nco/4.9.5

# MD July 8 2020
# Edit to form 2018 file only

set -x
set -e

if [ `ncdump -tv time work/2018123112.nc | grep 2019-01-01 | wc -l` -eq 0 ] ; then
  # last day is repeated, adjust it's time forward 1dh
  ncap2 -s "time=time+24" work/2018123112.nc work/tmp.nc && mv work/tmp.nc work/2018123112.nc
fi
 
# Concatenate into one long record
rm -f rdrs.nc
ncrcat -h work/201*.nc -o rdrs.nc

#Cut out 2018
ncks -F -d time,13,8772 rdrs.nc RDRSv21_y2018.nc
ncks -F -d time,8773,8784 rdrs.nc RDRSv21_y2019.nc

ncatted -a _FillValue,,o,f,NaN RDRSv21_y2018.nc
ncatted -a _FillValue,,o,f,-999.0 RDRSv21_y2018.nc
  
ncatted -a _FillValue,,o,f,NaN RDRSv21_y2019.nc
ncatted -a _FillValue,,o,f,-999.0 RDRSv21_y2019.nc
