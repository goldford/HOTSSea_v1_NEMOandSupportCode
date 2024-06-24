#!/bin/ksh
# in mobatexteditor ensure Unix format button is checked 
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;
#
# GO 2023-04-23
# clip ORAS data to JdF mouth
# this script meant to be run on unix (server: graham)
# downloaded 2024-03-05
# github repo does not include ORAS5 (854 mb) so skip to next script
yrs={1979..2018}
mos={"01","02","03","04","03","04","05","06","07","08","09","10","11","12"}

x_min=26
x_max=30
y_min=25
y_max=28

pathin="/project/6006412/mdunphy/Forcing/ORAS5/opa0_truncate/"
pathin_m="/project/6006412/mdunphy/Forcing/ORAS5/"
pathout="/project/6006412/goldford/data/atmos/ORAS5_JdF/"
pathtemp="/project/6006412/goldford/data_temp/ORAS5/"

dir="${pathtemp}" ;
if [[ ! -e $pathtemp ]] ; then
  mkdir $pathtemp
fi

echo clipping files to bounding box;
for y in $yrs ; do
    echo y: $y ;
    #pathin2="/CDF/"$y"/"
   # pathin3=$pathin1""$pathin2
  for m in ${mos[@]}; do
    echo ${pathin}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea.nc ;
    #votempervosaline_ORAS5_1m_197912_grid_T_02_salishsea.nc
    #cdo -L -sellonlatbox,${lon_min},${lon_max},${lat_min},${lat_max} -selname,votemper,vosaline ${pathin}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea.nc ${pathtemp}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea_JdF.nc 
    cdo -L -selindexbox,${x_min},${x_max},${y_min},${y_max} -selname,votemper,vosaline ${pathin}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea.nc ${pathtemp}votempervosaline_ORAS5_1m_${y}${m}_grid_T_02_salishsea_JdF.nc

  done
done

#cdo -L -selindexbox,${x_min},${x_max},${y_min},${y_max} ${pathin_m}mesh_mask_trunc.nc ${pathout}mesh_mask_trunc_JdF.nc
cdo -L -selindexbox,${x_min},${x_max},${y_min},${y_max} ${pathin_m}bathy_meter_trunc.nc ${pathout}bathy_meter_trunc_JdF.nc


echo concatinating files ;
dir="${pathout}" ;
if [[ ! -e $dir ]] ; then
  mkdir $dir
fi
ncrcat ${pathtemp}*ORAS*.nc ${pathout}votempervosaline_1m_ORAS5_1979_2018_JdF.nc
pathin="/project/6006412/goldford/data/atmos/ORAS5_JdF/"

echo interpolating to daily;
cdo -inttime,1979-01-16,12:00,1day ${pathout}votempervosaline_1m_ORAS5_1979_2018_JdF.nc ${pathout}votempervosaline_1d_ORAS5_1979_2018_JdF_interp.nc

