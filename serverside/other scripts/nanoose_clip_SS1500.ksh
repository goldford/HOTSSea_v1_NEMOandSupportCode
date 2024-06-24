#!/bin/ksh
# in mobatexteditor ensure Unix format button is checked 
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;
#
# GO 2023-0#3-22
# extract temp and salin from SS1500 to match Nanoose Stn obs (clip to bounding box using cdo)

yrs={1979..2018}
mos={"01","02","03","04","03","04","05","06","07","08","09","10","11","12"}
mod_run="RUN216"


# bounding box clip will result in  
# 6 x 5 rectangle (model grid indices: 71 <= i <= 76, 172 <= j <= 176).

# can also clip using indices
# cdo -L -selindexbox,${x_min},${x_max},${y_min},${y_max} -selname,votemper,vosaline ${pathin}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea.nc ${pathtemp}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea_JdF.nc

# this is rough bounding box from map
#lat_min=49.25 
#lat_max=49.449
#lon_min=-124.2
#lon_max=-123.749

# this is the centroid of pts (49.33, -124) plus minus a small amount
lat_min=49.31
lat_max=49.35
lon_min=-124.05
lon_max=-123.95

pathin1="/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-${mod_run}"
pathout="/project/6006412/goldford/data/evaluation/nanoose_stn/${mod_run}-SS1500-extract/"
pathtemp="/project/6006412/goldford/data_temp/clip_results/"

echo $pathtemp ;
dir="${pathtemp}" ;
if [[ ! -e $pathtemp ]] ; then
  mkdir $pathtemp
fi


echo clipping files to bounding box;

for y in $yrs ; do
    echo y: $y ;
    #pathin2="/CDF/"$y"/"
    pathin2="/CDF/"
    pathin3=$pathin1""$pathin2
  for m in ${mos[@]}; do
    echo ${pathin3}SalishSea1500-${mod_run}_1h_grid_T_y${y}m${m}.nc ;
    cdo -L -sellonlatbox,${lon_min},${lon_max},${lat_min},${lat_max} -selname,votemper,vosaline ${pathin3}SalishSea1500-${mod_run}_1h_grid_T_y${y}m${m}.nc ${pathtemp}SS1500-${mod_run}_1h_grid_T_y${y}m${m}_clipped.nc
  done
done


# interpolate to daily to reduce size - they will be averaged to biweekly anyway
echo concatinating files ;

ncrcat ${pathtemp}*SS1500-${mod_run}_1h_grid_T*.nc ${pathtemp}votempervosaline_SS1500-${mod_run}_1979_2018_nanstn.nc

echo interpolating to daily;
cdo -inttime,1979-01-16,12:00,1day ${pathtemp}votempervosaline_SS1500-${mod_run}_1979_2018_nanstn.nc ${pathtemp}votempervosaline_1d_${mod_run}_1979_2018_nanstn_interp.nc

echo $pathout ;
dir="${pathout}" ;

if [[ ! -e $pathout ]] ; then
  mkdir $pathout
fi
# compute average across all grid cells
echo computing avg over all grid cells;
ncwa -a y,x ${pathtemp}votempervosaline_1d_${mod_run}_1979_2018_nanstn_interp.nc ${pathout}votempervosaline_1d_${mod_run}_1979_2018_nanstn_avg.nc
