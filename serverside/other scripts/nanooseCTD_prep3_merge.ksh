#!/bin/ksh
#source $LMOD_PKG/init/ksh
#module load cdo/1.9.8 
#module load nco/4.9.5

# by G Oldford 2022-2023
# purpose: use cdo to merge nanoose stn data files (1 / cast) to annual files
#          then merge those together and calculate daily means across all years
# note: can't merge all at once - too many files open at once

pathin1="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_interp/"
pathout1="/project/6006412/goldford/data/evaluation/nanoose_stn/_temp/"
pathout2="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_mrgd/"
pathout3="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_mrgd_mean/"

echo $pathin1 ;
echo $pathout1 ;

if [[ ! -e $pathout1 ]] ; then
  mkdir $pathout1
fi

if [[ ! -e $pathout2 ]] ; then
  mkdir $pathout2
fi

if [[ ! -e $pathout3 ]] ; then
  mkdir $pathout3
fi

year_start=1981
year_end=1981

#for y in {${year_start}..${year_end}} ; do 
#  echo "pre-merging "${y} ;
#  file_in=${pathin1}*${y}*.nc
#  file_out=${pathout1}merged_${y}.nc
#
#  echo ${file_out} ;
#  
#  cdo -O mergetime ${file_in} ${file_out}
#  
#done

# merge indiv files to multi year
cdo -O mergetime ${pathout1}*{${year_start}..${year_end}}.nc ${pathout2}merged_${year_start}to${year_end}.nc

# daily mean across all years for climatology
cdo -O -ydaymean ${pathout2}merged_${year_start}to${year_end}.nc ${pathout3}ydaymean_${year_start}to${year_end}.nc
cdo -O -daymean ${pathout2}merged_${year_start}to${year_end}.nc ${pathout3}daymean_${year_start}to${year_end}.nc
