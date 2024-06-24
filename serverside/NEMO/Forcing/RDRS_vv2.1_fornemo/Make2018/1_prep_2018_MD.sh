#!/bin/bash
module load nco/4.9.5

# MD July 8 2020
# Edit to form 2018 file only

set -x
set -e

mkdir -p work

for fin in /home/mdunphy/project/goldford/data_temp/RDRSv21/rdrs2000son/{20171231,2018}*.nc; do
  echo $fin

  rm -f temp.nc

  # Rename to less awkward names 
  ncrename -v RDRS_v2.1_P_UUC_10m,u10m \
           -v RDRS_v2.1_P_VVC_10m,v10m \
           -v RDRS_v2.1_P_FI_SFC,therm_rad \
           -v RDRS_v2.1_P_FB_SFC,solar_rad \
           -v RDRS_v2.1_A_PR0_SFC,precip \
           -v RDRS_v2.1_P_PR0_SFC,precip_fcst \
           -v RDRS_v2.1_P_PN_SFC,slp \
           -v RDRS_v2.1_P_HU_1.5m,h1p5m \
           -v RDRS_v2.1_P_TT_1.5m,t1p5m \
           $fin temp.nc

  # Drop unneeded variables
  ncks -C -x -v RDRS_v2.1_P_GZ_SFC,RDRS_v2.1_P_TD_1.5m,RDRS_v2.1_P_TD_09944,RDRS_v2.1_P_HU_09944,RDRS_v2.1_P_TT_09944,RDRS_v2.1_P_GZ_09944,RDRS_v2.1_P_GZ_SFC,rotated_pole,rlat,rlon temp.nc -O temp.nc

  # Change units and units attribute
  ncap2 -s u10m=u10m*0.514444 \
        -s v10m=v10m*0.514444 \
        -s t1p5m=t1p5m+273.15 \
        -s slp=slp*100        \
        -s precip=precip*1000/3600 \
        temp.nc -O -o temp.nc

  # Compute snow
  ncap2 -s "snow=float(precip);where(t1p5m>=273.15)snow=0.0" temp.nc -O -o temp.nc

  ncatted -a units,u10m,m,c,"m s**-1" \
          -a units,v10m,m,c,"m s**-1" \
          -a units,t1p5m,m,c,"Kelvin" \
          -a units,slp,m,c,"Pa"       \
          -a units,precip,m,c,"kg m**-2 s**-1" \
          -a units,snow,m,c,"kg m**-2 s**-1" \
          -O temp.nc

  fout=`basename $fin`
  mv temp.nc work/$fout

done

