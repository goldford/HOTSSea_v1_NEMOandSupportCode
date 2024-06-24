#!/bin/ksh

set -ex

## ======================================================================
## AUTO-RESTART: PRE-PROCESSING SCRIPT
## RAPPAT.ksh
## ======================================================================

# Import run settings
. ./environment.ksh

## ======================================================================
## PREPARE THE RUN 
## ======================================================================

if [ -z $DIR_TMP_RUN ] ; then
  echo "DIR_TMP_RUN IS NOT DEFINED, CAN NOT PROCEED"
  exit
fi
cd $DIR_TMP_RUN
rm -rf $DIR_TMP_RUN/*	# Clean DIR_TMP_RUN in automated sequence

## ======================================================================
## Copy the binaries
## ======================================================================

if [ ! -f ${DIR_CALCU_EXE}/$EXEC_NAME ] ; then
    echo ${DIR_CALCU_EXE}/$EXEC_NAME 'IS MISSING !!!!!'
    exit
else
    cp ${DIR_CALCU_EXE}/$EXEC_NAME ./opa
    chmod 755 ./opa
fi
# --- XIOS SERVER
if [ $NBPROC_XIO -gt 0 ] ; then
    cp ${DIR_CALCU_XIO}/xios_server.exe ./xios_server.exe
    chmod 755 ./xios_server.exe
fi

## ======================================================================
## Copy the config xml files 
## ======================================================================
if [ ! -f ${DIR_CALCU_CTL}/iodef-${CONFIG}.xml ] ; then
    echo ${DIR_CALCU_CTL}/iodef-${CONFIG}.xml 'IS MISSING !!!!!'  ;  exit
else
    cp ${DIR_CALCU_CTL}/iodef-${CONFIG}.xml iodef.xml
fi

if [ ! -f ${DIR_CALCU_CTL}/domain_def-${CONFIG}.xml ] ; then
    echo ${DIR_CALCU_CTL}/domain_def-${CONFIG}.xml 'IS MISSING !!!!!' ; exit
else
    ln -s ${DIR_CALCU_CTL}/domain_def-${CONFIG}.xml domain_def.xml
fi

if [ ! -f ${DIR_CALCU_CTL}/file_def-${CONFIG}.xml ] ; then
    echo ${DIR_CALCU_CTL}/file_def-${CONFIG}.xml ' is missing, we must be using xios-1, skipping'
else
    cp    ${DIR_CALCU_CTL}/file_def-${CONFIG}.xml file_def.xml
fi

if [ ! -f ${DIR_CALCU_SHR}/field_def.xml ] ; then
    echo ${DIR_CALCU_SHR}/field_def.xml 'IS MISSING !!!!!' ;  exit
else
    ln -s ${DIR_CALCU_SHR}/field_def.xml .
fi

# ======================================================================
## Copy the namelists
## ======================================================================
if [ ! -f $DIR_CALCU_CTL/namelist_cfg-${CONFIG} ] ; then
    echo $DIR_CALCU_CTL/namelist_cfg-${CONFIG} 'IS MISSING !!!!!' ;   exit
else
    cp $DIR_CALCU_CTL/namelist_cfg-${CONFIG} namelist_cfg
fi

if [ ! -f $DIR_CALCU_SHR/namelist_ref ] ; then
    echo $DIR_CALCU_SHR/namelist_ref 'IS MISSING !!!!!' ;    exit
else
    cp $DIR_CALCU_SHR/namelist_ref namelist_ref
fi

## ======================================================================
## Determine the nit000, nitend for this run and the start date
## ======================================================================

# Add up total number of time steps completed so far
niter=0
while read no nit000 nitend; do
  niterc=$niter               # number of time steps completed
  dif=$(($nitend - $nit000 + 1))
  niter=$(($niter + $dif))    # number of time steps described in file
done < $DIR_CALCU_CTL/$CONFIG.db

# Reload last line of database
read no nit000 nitend < <(tail -1 $DIR_CALCU_CTL/$CONFIG.db)

# Determine start date for this run
if [ $no -eq 1 ]; then
  ndate0=$NDATEDEB
  clstep=0
  if [ $RESTART == "yes" ]; then
      flag_rst=.true.
      flag_ramp=.false.
  else
      flag_rst=.false.
      flag_ramp=.true.
  fi
  rsttype=0
else
  let "days_completed = $niterc * $PASDETEMPS / 86400"
  ndate0=$(/bin/date "+%Y%m%d" --date="${NDATEDEB} ${days_completed} days")
  # Find the time step number where the previous run ended
  read discard discard clstep < <(tail -2 $DIR_CALCU_CTL/$CONFIG.db | head -1)
  flag_rst=.true.
  flag_ramp=.false.
  rsttype=2    
fi

## ======================================================================
## Check consistency between SPLIT_FREQ (d), NNSTOCK (ts) & total 
## number of time steps for the time slice 
## ======================================================================
checksplit=$((${NNSTOCK}/(86400/${PASDETEMPS}) ))   
totTS=$(( $nitend - $nit000 + 1 ))
TSperday=$(( 86400/${PASDETEMPS}  ))
goodtogo=1

if [ $checksplit -ne $SPLIT_FREQ ] ; then
  echo ' W A R N I N G : checksplit NOT EQUAL TO SPLIT_FREQ changes were done to environment.ksh' 
  goodtogo=0
fi

modts=$(( $totTS % $TSperday ))
if [ $modts -ne 0 ] ; then
  echo ' W A R N I N G : total number of timestep not a factor of # timestep per day'
  goodtogo=0
fi

stockmod=$(( ${NNSTOCK} % $TSperday ))
if [ $stockmod -ne 0 ] ; then 
  echo ' W A R N I N G : NNSTOCK not a factor of # timestep per day'
  goodtodo=0
fi

modstock=$(( $totTS % $NNSTOCK ))
if [ $modstock -ne 0 ] ; then 
  echo ' W A R N I N G : totTS not a factor of NNSTOCK'
  goodtogo=0
fi

tspersplit=$(( $SPLIT_FREQ * (86400/${PASDETEMPS})  ))
checksplit=$(( $totTS % tspersplit ))
if [ $modstock -ne 0 ] ; then        
  echo ' W A R N I N G : split_freq not a factor of totTS'
  goodtogo=0
fi

if [ ${MESHRUN} = 'yes' ] ; then 
  echo ' W A R N I N G  RUN FOR MESH MASK !!'
  echo ' CONSISTENCY CHECK OVERWRITTEN'
  goodtogo=1
fi

if [ $goodtogo == 0 ] ; then
  echo ' E R R O R : check your consistency SPLIT_FREQ (d), NNSTOCK (ts)'
  echo '             & total number of time steps for the time slice '
  echo ' A B O R T !'
  exit
else 
  echo ' G O O D   T O   G O !!!'
fi


## ======================================================================
## Calendar
## ======================================================================
ndays=` echo 1 | awk "{ a=int( ($nitend - $nit000 +1)*$PASDETEMPS /86400.) ; print a }"`
ydeb=`date --date="${ndate0}" +%Y`
mdeb=`date --date="${ndate0}" +%m`
ddeb=`date --date="${ndate0}" +%d `
ndastpfin=$(/bin/date "+%Y%m%d" --date="${ndate0} ${ndays} days")
yfin=`date --date="${ndastpfin}" +%Y`
mfin=`date --date="${ndastpfin}" +%m`
dfin=`date --date="${ndastpfin}" +%d `
julstart=$(d2j ${ndate0})
julstop=$(d2j ${ndastpfin})
echo $ndays days to run, starting $ndate0 "("$julstart")" ending $ndastpfin "("$julstop")"

## ======================================================================
## Rappat : Copy the files adapted to the config
## ======================================================================
## ======================================================================
## Bathymetry
## ======================================================================
nomini=$BATMET_FIL
nom=bathy_meter.nc
if [ ! -z ${nomini} ] ; then
   ln -sf ${DIR_ARCHI_INPUT}/${nomini} ${DIR_TMP_RUN}/${nom}
else
   echo ${DIR_ARCHI_INPUT}/${nomini} ' not found' ; exit
fi
## ======================================================================
## Coordinate
## ======================================================================
nomini=$COORD_FIL
nom=coordinates.nc
if [ ! -z ${nomini} ] ; then
   ln -sf ${DIR_ARCHI_INPUT}/${nomini} ${DIR_TMP_RUN}/${nom}
else 
   echo ${DIR_ARCHI_INPUT}/${nomini} ' not found' ; exit
fi
## ======================================================================
## Bottom friction 
## ======================================================================
#nomini=$BFRI_FIL
#nom=bfr_coef.nc
#if [ ! -z ${nomini} ] ; then
#   ln -sf ${DIR_ARCHI_INPUT}/${nomini} ${DIR_TMP_RUN}/${nom}
#fi
## ======================================================================
## Runoff
## ======================================================================
nomini='runoff_daigreen' #$RUNOFF_FIL
if [ $nomini != 'runoff_daigreen' ] ; then       # CLIMATOLOGYCAL RUNOFFS
  nom='runoff_1m_nomask.nc'
  ln -sf $DIR_ARCHI_INPUT/${nomini} ${DIR_TMP_RUN}/${nom}
fi

if [ ! -s $RUNOFF_FIL -o ! -z $RUNOFF_TEMP_FIL -o ! -z $RUNOFF_MASK_FIL ] ; then 
  mkdir -m 755 -p ${DIR_TMP_RUN}/RIVER
fi

if [ ! -z ${RUNOFF_FIL} ] ; then
    ln -sf ${DIR_ARCHI_INPUT}/${RUNOFF_FIL}y${yy}*.nc       ${DIR_TMP_RUN}/RIVER/.
fi
if [ ! -z ${RUNOFF_TEMP_FIL} ] ; then
    ln -sf ${DIR_ARCHI_INPUT}/${RUNOFF_TEMP_FIL}y${yy}*.nc  ${DIR_TMP_RUN}/RIVER/.
fi
if [ ! -z ${RUNOFF_MASK_FIL} ] ; then
    ln -sf ${DIR_ARCHI_INPUT}/${RUNOFF_MASK_FIL}.nc         ${DIR_TMP_RUN}/RIVER/.
fi


## ======================================================================
## physics obc files 
## ======================================================================
if [ $OBC_ACTIVE = 'yes' ] ; then
    if [ ! -d $DIR_TMP_OBC ] ; then
       mkdir -p $DIR_TMP_OBC
    fi
    for obc in ${obc_list} ; do
      obcin=$DIR_ARCHI_OBC/$obcpref2D'_'$obcsuf2D'.nc'
      if [ ! -f $obcin ] ; then
         echo $obcin 'IS MISSING !!!!!' 
         exit
      else
         ln -fs $obcin $DIR_TMP_OBC/.
      fi

      if [ ! -z $obcprefNorth ] ; then 
         obcN=$DIR_ARCHI_OBCN/$obcprefNorth'_*.nc'
         if [ ! -f $obcN ] ; then
            echo $obcN 'IS MISSING !!!!!' 
            exit
         else
            ln -fs $obcN $DIR_TMP_OBC/.
         fi
      fi
    done
fi

## ==================================================================
## physics tidal files - Added : FCH
## ==================================================================
if [ $TIDE_ACTIVE = 'yes' ]; then
  if [ $TIDE_2Dfile = 'yes' ]; then
     if [ ! -s $DIR_ARCHI_TIDE/${TIDEH_FIL} -o \
          ! -s $DIR_ARCHI_TIDE/${TIDEU_FIL} -o \
          ! -s $DIR_ARCHI_TIDE/${TIDEV_FIL} ] ; then 
        echo ' SOME 2D TIDE FILES MISSING ' ; exit
     fi
     # tide name hardcoded with prefix defined in namelist.
     ln -sf $DIR_ARCHI_TIDE/${TIDEH_FIL} ${DIR_TMP_RUN}/tide_grid_T.nc
     ln -sf $DIR_ARCHI_TIDE/${TIDEU_FIL} ${DIR_TMP_RUN}/tide_grid_U.nc
     ln -sf $DIR_ARCHI_TIDE/${TIDEV_FIL} ${DIR_TMP_RUN}/tide_grid_V.nc

   else
     echo ' YOU SHOULD NOT GO THERE ' ; exit
#
#     grid_list='T U V'
#     for obc in ${obc_list} ; do
#       for compo in ${compo_list}; do
#         for grid in ${grid_list}; do
#            #nom=$TIDE_PREF'_'$obc'_'$compo'_grid_'$grid'.nc'#LUO
#            # nomini=$TIDE_PREF'_'$obc'_grid_'$grid'.nc'  #luo
#            nom=$TIDE_PREF'_grid_'$grid'.nc' #luo
#            # nomini=$TIDE_PREF'_grid_'$grid'.nc' #luo
#            tidein=$DIR_ARCHI_TIDE/$nom #LUO
#            if [ ! -f $tidein ] ; then
#              echo $tidein 'IS MISSING !!!!!' 
#              exit
#            else
#               ln -fs $tidein $DIR_TMP_OBC/$nom
#             fi
#         done
#       done
#     done
   fi
   #Zhao Peng tide ramp is only used at cold start
    TIDERAMP=.true.
    if [ $no -ne 1 ]; then 
       TIDERAMP=.false.  
       echo $TIDERAMP $flag_rst 'flagrst'
    fi
  echo $TIDERAMP $flag_rst 'tideramp'
else
   TIDERAMP=.false. # default
fi

#### B02 ####
#TIDERAMP=.false.
#flag_ramp=.false.

## ======================================================================
## Restart
## ======================================================================

if [ ! -d ${DIR_TMP_RUN}/OUTRST ] ; then mkdir -m 755 -p ${DIR_TMP_RUN}/OUTRST ; fi

if [ "$RESTART" == "no" ] && [ $no == 1 ] ; then
     #ln -s $DIR_ARCHI_CLIM/$TEMP_FIL $DIR_TMP_RUN/IC_T.nc   # Temperature   
     #ln -s $DIR_ARCHI_CLIM/$SAL_FIL  $DIR_TMP_RUN/IC_S.nc   # Salinity
     #ln -s $DIR_ARCHI_CLIM/$U_FIL    $DIR_TMP_RUN/IC_U.nc   # U vel
     #ln -s $DIR_ARCHI_CLIM/$V_FIL    $DIR_TMP_RUN/IC_V.nc   # V vel
     #ln -s $DIR_ARCHI_CLIM/$SSH_FIL  $DIR_TMP_RUN/IC_SSH.nc # SSH
    echo "This is a spinup run"
else
    typeset -Z8 clstep
    nomrest=restart.nc
    nomrestart=${CONFIG}_${clstep}_restart.nc
    if [ $no -eq 1 ] && [ -s $RESTART_FILE ] ; then
        # Look for combined restart file with file name specified in $RESTART_FILE
        echo "FIRST RUN, COMBINED RESTART FILE FOUND: " $RESTART_FILE
        ln -s $RESTART_FILE $nomrest
    elif [ -s $DIR_ARCHI_RST/$nomrestart ] ; then
        # Look for combined restart file with expected file name
        echo "FOUND COMBINED RESTART FILE: " $DIR_ARCHI_RST/$nomrestart
        ln -s $DIR_ARCHI_RST/$nomrestart $nomrest
    else
        echo "DID NOT FIND COMBINED RESTART FILE, LOOKING FOR PER-PROCESSOR RESTART FILES"
        # JPP : Restart file names (or prefix) are defined in
        #       namelist (cn_ocerst_in) and namelist_ice (cn_icerst_in)
        #       and must fit the prefix here!
        #       Nb of Proc from 0 to NBPROC-1, jproc must be 4 digits
        typeset -Z4 jproc
        for jproc in `seq 0 $(($NBPROC-1))`; do
            nomrest=restart_${jproc}.nc
            nomrestart=${CONFIG}_${clstep}_restart_${jproc}.nc
            if [ ! -s $DIR_ARCHI_RST/$nomrestart ] ; then
                echo ' MISSING RESTART FILE ' $DIR_ARCHI_RST/$nomrestart ; exit
            fi
            ln -s $DIR_ARCHI_RST/$nomrestart $nomrest

#            # SEA ICE -- JPP NOT CHECKED NEMO3.6 (no CICE in our version)
#            if [ $SEAICE_MODULE == 'LIM' || $SEAICE_MODULE == 'CICE' ]; then
#                nomrest=restart_ice_${jproc}.nc
#                nomrestart=${CONFIG}_${clstep}_restart_ice_${jproc}.nc
#                ln -s $DIR_ARCHI_RST/$nomrestart $nomrest
#            fi
        done
    fi
fi

## ======================================================================
## Atmospheric forcings
## ======================================================================
if [ $ATMOS_DATA != 'NONE' ] ; then

  if [ ! -d $DIR_TMP_ATM ] ; then mkdir $DIR_TMP_ATM ; fi
  
  # Start links at previous day
  datec=$(/bin/date "+%Y%m%d" --date="${ndate0} -1 days")

  # Flag to proceed if first day is missing (can occur when we start on the first day of atm forcing)
  fmok=1

  while [ $datec -le $ndastpfin ] ; do 
    year=`date  --date="${datec}" +%Y`
    month=`date --date="${datec}" +%m`
    day=`date   --date="${datec}" +%d `

    if [ $ATMOS_DATA == 'HRDPS-WEST' ] ; then
      prefix='ops'
      nomfic=${prefix}_y${year}m${month}d${day}.nc
    elif [ $ATMOS_DATA == 'HRDPS-EAST' ] ; then
      prefix='meopar_east'
      nomfic=${year}/${prefix}_${year}${month}${day}.nc
    elif [ $ATMOS_DATA == 'OPPwest' ] ; then
      prefix='HRDPS_OPPwest_ps2.5km'
      nomfic=${year}/${prefix}_y${year}m${month}d${day}.nc
    elif [ $ATMOS_DATA == 'OPPwestSS' ] ; then
      prefix='HRDPS_OPPwestSS_ps2.5km'
      nomfic=${year}/${prefix}_y${year}m${month}d${day}.nc
    fi

    nomfile=${prefix}_y${year}m${month}d${day}.nc

    if [ $fmok == 1 ] && [ ! -f $DIR_ARCHI_ATMOS/${nomfic} ] ; then 
      echo "Previous day atmos file $nomfile not found, proceeding anyway"
    elif [ ! -f $DIR_ARCHI_ATMOS/${nomfic} ] ; then
      echo $DIR_ARCHI_ATMOS/${nomfic} ' DOES NOT EXIST ... EXIT' ; exit
    else
      ln -s $DIR_ARCHI_ATMOS/${nomfic} ${DIR_TMP_ATM}/${nomfile}
    fi

    # Advance date by one day for next trip through loop
    datec=$(/bin/date "+%Y%m%d" --date="${datec} 1 days")

    # Do not permit any more days to be missing
    fmok=0
  done

  nomweight='weights.nc'
  if [ ! -f $DIR_ARCHI_INPUT/$WEIGHTS_FIL ] ; then 
    echo $DIR_ARCHI_INPUT/$WEIGHTS_FIL ' DOES NOT EXIST ... EXIT' ; exit
  else
    ln -s $DIR_ARCHI_INPUT/$WEIGHTS_FIL ${DIR_TMP_ATM}/${nomweight}
  fi
fi

## ======================================================================
## Salish Sea forcing
## ======================================================================
if [ "$CONF_NAME" == "SalishSea1500" ] ; then
  # BDY from ORAS5 W, Dosser N
  src=/home/mdunphy/project/mdunphy/Forcing/forcing_oras5_20230717   # third attempt at bias fix
  ln -s ${src}/bdy_ssh bdy_ssh
  ln -s ${src}/bdy_ts bdy_ts
  ln -s ${src}/bdy_uv bdy_uv
  ln -s /home/mdunphy/project/mdunphy/Forcing/Dosser_north_bdy/Dosser_SalishSea1500_N.nc
  
  # IC from S. Allen
  src=/home/mdunphy/project/mdunphy/Forcing/sallen_ts
  ln -s ${src}/ic_ts ic_ts

  # Atmos forcing RDRS v2.1 prepped by G. Oldford
  ln -s /home/mdunphy/project/mdunphy/Forcing/RDRS_v2.1_fornemo RDRSv21

  ln -s /home/mdunphy/project/mdunphy/Forcing/ss1500tides tides
  ln -s /project/6049289/goldford/data/runoff/rivers_month_20210406_GO.nc .
#  ln -s /home/mdunphy/project/mdunphy/Forcing/daily_runoff_1979_2019 rivers  # path changed 
  ln -s /home/mdunphy/project/mdunphy/Forcing/daily_runoff_1979_2019_newpath rivers  # path changed 
  ln -s /home/mdunphy/project/mdunphy/Forcing/grid/no_ice.nc .
fi

## ======================================================================
## namelist updates
## ======================================================================
sed -i \
    -e "s/NUMERO_DE_RUN/$no/"        \
    -e "s/EXPER/\"$CONFIG\"/"        \
    -e "s/NIT000/$nit000/"           \
    -e "s/NITEND/$nitend/"           \
    -e "s/NDATE0/$ndate0/"           \
    -e "s/BREST/$flag_rst/"          \
    -e "s/NNSTOCK/$NNSTOCK/"         \
    -e "s/JPNIJ/$NBPROC/"            \
    -e "s/JPNI/$JPNI/"               \
    -e "s/JPNJ/$JPNJ/"               \
    -e "s/TIDERAMP/$flag_ramp/"      \
    -e "s/DATERST/$rsttype/"         \
    -e "s/PASDETEMPS/$PASDETEMPS/"   \
    -e "s/NNBARO/$NNBARO/"           namelist_cfg

## Final modification to the iodef.xml -- to includde the split_freq in days
sed -i -e "s/SPLIT_FREQ/${SPLIT_FREQ}d/" iodef.xml


## ======================================================================
## Soumission
## ======================================================================
cd $DIR_CALCU_CTL

echo "Submitting job, ${TOTPROC} cores, ${TIMELIMIT}s runlimit, control path $DIR_CALCU_CTL"
#exit
#ord_soumet $CONF_NBR.ksh  -image ${CONTAINER} -mach ${MACHINE} -cpus ${TOTPROC} -mpi -cm ${MEMLIMIT} \
#                          -t ${TIMELIMIT} -listing ${DIR_CALCU_CTL} -q ${QUEUE} \
#                          -o ${DIR_CALCU_CTL} -project ${PROJECT} \
#                          -mail Michael.Dunphy@dfo-mpo.gc.ca -notify -keep
ACCOUNT=''
#ACCOUNT=--account=${PROJECT}
#echo sbatch --time=${TIMELIMIT} --job-name=$CONF_NBR --mem-per-cpu=${MEMLIMIT} --ntasks=${TOTPROC} ${ACCOUNT} -o ${DIR_TMP_RUN}/nemo.stdout -e ${DIR_TMP_RUN}/nemo.stderr $CONF_NBR.ksh
#     sbatch --time=${TIMELIMIT} --job-name=$CONF_NBR --mem-per-cpu=${MEMLIMIT} --ntasks=${TOTPROC} ${ACCOUNT} -o ${DIR_TMP_RUN}/nemo.stdout -e ${DIR_TMP_RUN}/nemo.stderr $CONF_NBR.ksh
echo sbatch --time=${TIMELIMIT} --job-name=$CONF_NBR --mem-per-cpu=${MEMLIMIT} -N 3 --ntasks=${TOTPROC} ${ACCOUNT} -o ${DIR_TMP_RUN}/nemo.stdout -e ${DIR_TMP_RUN}/nemo.stderr $CONF_NBR.ksh
sbatch --time=${TIMELIMIT} --job-name=$CONF_NBR --mem-per-cpu=${MEMLIMIT} -N 3 --ntasks=${TOTPROC} ${ACCOUNT} -o ${DIR_TMP_RUN}/nemo.stdout -e ${DIR_TMP_RUN}/nemo.stderr $CONF_NBR.ksh

## ======================================================================
## End of the rappat script
## ======================================================================
