#!/bin/ksh
#source $LMOD_PKG/init/ksh
#module load cdo/1.9.8 
#module load nco/4.9.5
# by G Oldford 2022-2023
# processes separate prepped Nanoose std NC files to a single file.
# command prompt example: 'ksh prep_nstnCTD_merge.ksh'

pathin1="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped/"
echo $pathin1 ;

pathtemp="/project/6006412/goldford/data_temp"
pathout="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_singlefile"
pathout1="/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_singlefile/"

echo $pathout ;
dir="${pathout}" ;

if [[ ! -e $dir ]] ; then
  mkdir $dir
fi

# calc monthly mean of already daily averaged fields stored in monthly NC files
echo concatinating files ;

ncrcat ${pathin1}*.nc ${pathout1}Nanoose_alltgthr.nc