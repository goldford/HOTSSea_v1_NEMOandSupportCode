#!/bin/ksh
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;

# GO 2022-02-12
# lighthouse data from multiple years needs to get mashed together

pathin="/project/6006412/goldford/data/evaluation/Lighthouses/processed/annual_files/"
pathout="/project/6006412/goldford/data/evaluation/Lighthouses/processed/"
#pathtemp="../DATA/TEMP/" # for intermed files

year_start=1979
year_end=2019

# Declare an array of string with type
declare -a StringArray=("active_pass_lightstation" "amphitrite_point_lightstation" "cape_beale_lightstation" 
"cape_mudge_lightstation" "chrome_island_lightstation" "departure_bay_(pbs)" "entrance_island_lightstation" 
"race_rocks_lightstation" "sheringham_point_lightstation" "sisters_islets_lightstation" "west_van_labs_(caer)" 
"bamfield_marine_sciences_centre")

#declare -a StringArray=("entrance_island_lightstation" )
 
# Iterate the string array using for loop
for val in ${StringArray[@]}; do
  echo ${pathin}${val}_SalinTemp_;
  ncrcat ${pathin}${val}_SalinTemp_* ${pathout}${val}_SalinTemp_${year_start}_${year_end}.nc
done


#for y in {${year_start}..${year_end}} ; do 

#  file="file$i"
#  dest="dest$i"
#  echo "${!file} ${!dest}" #or whatever you want to do with each file
#done
#done