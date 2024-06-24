#!/bin/ksh
set -ax
##	    environment.ksh

## ======================================================================
## ===============
## 2014/03/05 : JP Paquin
## Defines the environment and PATH for running NEMO with automatic 
## resubmission capacities for NEMO 3.1
## HEAVILY BASED ON MERCATOR VERSION
## ===============
## Revision: Fatemeh - Adapted for NEMO 3.6
##           JP Paquin - 2017/09/03 - refit and clean based on OPP BoF180...
##                                    to be used in OPP first and extended later.
## Revision: M.Dunphy, 2019/03, tidying/simplifying for OPP FA12
## ======================================================================

## =====================
## Execution environment
## =====================
## MACHINE       : cluster name (gpsc1, gpsc2 or gpsc4)
## CONTAINER     : container name (ex: dfo/dfo_all_default_ubuntu-18.04-amd64_latest)
## BACKEND       : script containing SSM pacages and env vars to use
## PROJECT       : GPSC project tag
## TIMELIMIT     : Time limit for NEMO runs, 14400=4h, 21600=6h
## TIMELIMIT_PP  : Time limit for port processing, 3600=1h
## MEMLIMIT      : Memory limit (per core, not per node)
## QUEUE         : Name of queue
## JPNI          : MPI decomposition along i (namelist JPNI)
## JPNJ          : MPI decomposition along j (namelist JPNJ)
## NBPROC        : number of NEMO processors (namelist JPNIJ)
## NBPROC_XIO    : number of XIOS processors in detached mode (0 if attached mode)
MACHINE=gpsc4
CONTAINER=dfo/dfo_all_default_ubuntu-18.04-amd64_latest
#BACKEND=backend_18.04_intel2019_static.ksh
PROJECT=rrg-nereusvc
TIMELIMIT=01:00:00
TIMELIMIT_PP=01:00:00
MEMLIMIT=3000M
QUEUE=dev

JPNI=6
JPNJ=23
NBPROC=94
NBPROC_XIO=2

#JPNI=22
#JPNJ=31
#NBPROC=338
#NBPROC_XIO=14
TOTPROC=`expr $NBPROC + $NBPROC_XIO`

## ==================
## Run length control
## ==================
## MAXSUB        : automatic resubmit up to MAXSUB lines in $CONFIG.db
## NDATEDEB      : date beginning of simulation
## PASDETEMPS    : Time step (in seconds); 86400/PASDETEMPS must be integer
## SPLIT_FREQ    : Frequency of writing output files (in days) for .xml file (for simplicity make it identical to NNSTOCK)
##                 -> code a check that it is a factor of the number of timesteps in the .db file in RAPPAT.ksh
## NNBARO        : NEMO3.6 ; nnbaro (# of baro sub-steps) now explicitely defined (NNBARO in namelist)
## NNSTOCK       : Restart writing frequency (default in restarting runs=99999999; otherwise select a factor of the number of days!)
## MESHRUN       : RUN just to get mesh mask; => Exception to verification of consistency
## RESTART       : 'yes' or 'no' -> Begin from restart files (if first run or not N.B. Harcoded names RAPPAT.ksh)
MAXSUB=26
NDATEDEB=20151129
PASDETEMPS=120
NNBARO=20
NNSTOCK=7200
#PASDETEMPS=20
#NNSTOCK=4320
SPLIT_FREQ=$((${NNSTOCK}/(86400/${PASDETEMPS}) ))               # Frequency of creation of output files
MESHRUN='no'
RESTART='no'

## ===============
## Run settings
## ===============
## CONF_NAME     : Script used to run the simulation
## CONF_NBR      :
## SEAICE_MODULE : Sea ice module to use (LIM/CICE/NONE)
## ATMOS_DATA    : Atmospheric forcing (HRDPS-EAST or HRDPS-WEST)
## OBC_ACTIVE    : Active OBC in the config
## TIDE_ACTIVE   : links to tidal files ()
## ======================================================================
CONF_NAME=SalishSea1500
CONF_NBR=RUN100
CONFIG=${CONF_NAME}-${CONF_NBR}
SEAICE_MODULE='NONE'
ATMOS_DATA='OPPwest'
OBC_ACTIVE='no'
TIDE_ACTIVE='no' # in SJAP100 the tidal signal is from the BDYs...
TIDE_2Dfile='no'
RIVER_ACTIVE='no'    # jpp-- sert a rien dans le code (?) 

## ======================================================================
## CALCU directories
## ======================================================================
## DIR_CALCU_EXE    : directory containing the opa binary         (rappat)
## DIR_CALCU_SHR    : directory containing the shared files       (rappat)
## DIR_CALCU_RBL    : directory containing the rebuild_nemo tool  (rappat)
## DIR_CALCU_XIO    : directory containing the xios server binary (rappat)
## DIR_CALCU_CTL    : directory containing control files (nml, xml, ..)
## DIR_CALCU_UTL    : directory where to find utilities
## ======================================================================
DIR_CALCU_EXE=${HOME}/projects/rrg-nereusvc/goldford/Software/NEMO36OPP/NEMOGCM/CONFIG/SalishSea1500/BLD/bin
DIR_CALCU_SHR=${HOME}/projects/rrg-nereusvc/goldford/Software/NEMO36OPP/NEMOGCM/CONFIG/SHARED
DIR_CALCU_RBL=${HOME}/projects/rrg-nereusvc/goldford/Software/NEMO36OPP/NEMOGCM/TOOLS/REBUILD_NEMO
DIR_CALCU_XIO=${HOME}/projects/rrg-nereusvc/goldford/Software/XIOS/bin
DIR_CALCU_CTL=${HOME}/projects/rrg-nereusvc/goldford/control/${CONF_NAME}-${CONF_NBR}
DIR_CALCU_UTL=${HOME}/bin # JPP not used at the moment (reserve for post-processing)

## ======================================================================
## ARCHIVE directories
## ======================================================================
## DIR_ARCHI_FORC      : directory to stock the forcing fields		(rappat)
## DIR_ARCHI_INPUT     : directory for the input files			(rappat)
## DIR_ARCHI_CLIM      : directory for the climatology files or IC      (rappat)
## DIR_ARCHI_OBC       : directory for the obc files			(rappat)
## DIR_ARCHI_TIDE      : directory for the tide files                   (rappat)
## DIR_ARCHI_OUTPUT    : directory to stock the results			(post_?)
## DIR_ARCHI_RST       : directory to stock the restarts		(post_?)
## ======================================================================
# Input files
DIR_ARCHI_FORC=${HOME}/projects/rrg-nereusvc/goldford/Forcing/grid
DIR_ARCHI_ATMOS=${HOME}/projects/rrg-nereusvc/mdunphy/Forcing/OPPWest
DIR_ARCHI_INPUT=${HOME}/projects/rrg-nereusvc/goldford/Forcing/grid
DIR_ARCHI_CLIM=/unused
DIR_ARCHI_OBC=/unused
DIR_ARCHI_OBCN=/unused
DIR_ARCHI_TIDE=/unused

# NEMO-generated files
DIR_ARCHI_HOME=${HOME}/projects/rrg-nereusvc/goldford/nemo_results/${CONF_NAME}/${CONFIG}
DIR_ARCHI_CDF=${DIR_ARCHI_HOME}/CDF
DIR_ARCHI_LOG=${DIR_ARCHI_HOME}/LOG
DIR_ARCHI_RST=${DIR_ARCHI_HOME}/RST

if [ ! -d $DIR_ARCHI_HOME ] ; then mkdir -m 755 -p $DIR_ARCHI_HOME ; fi
if [ ! -d $DIR_ARCHI_RST  ] ; then mkdir $DIR_ARCHI_RST ; fi
if [ ! -d $DIR_ARCHI_CDF  ] ; then mkdir $DIR_ARCHI_CDF ; fi
if [ ! -d $DIR_ARCHI_LOG  ] ; then mkdir $DIR_ARCHI_LOG ; fi


## ======================================================================
## TEMP directories
## =======================================================================
## DIR_TMP_HOME   : home temp directory 
## DIR_TMP_RST    : directory to put the restarts
## DIR_TMP_CDF    : directory to put the netcdfs
## DIR_TMP_LOG    : directory to put the logs
## DIR_TMP_RUN    : directory for run
## DIR_TMP_ATM    : directory for atmos symlinks
## DIR_TMP_OBC    : directory for obc fsymlinks
## DIR_TMP_RVR    : directory for river fsymlinks

## ======================================================================
DIR_TMP_HOME=${HOME}/projects/rrg-nereusvc/goldford/nemo_temp/${CONF_NAME}/${CONFIG}
DIR_TMP_RST=${DIR_TMP_HOME}/RST
DIR_TMP_CDF=${DIR_TMP_HOME}/CDF
DIR_TMP_LOG=${DIR_TMP_HOME}/LOG

DIR_TMP_RUN=${DIR_TMP_HOME}/RUN
DIR_TMP_ATM=${DIR_TMP_RUN}/ATMDATA
DIR_TMP_OBC=${DIR_TMP_RUN}/BDY
DIR_TMP_RVR=${DIR_TMP_RUN}/RIVER

if [ ! -d $DIR_TMP_HOME ] ; then mkdir -m 755 -p $DIR_TMP_HOME ; fi
if [ ! -d $DIR_TMP_RST  ] ; then mkdir $DIR_TMP_RST ; fi
if [ ! -d $DIR_TMP_CDF  ] ; then mkdir $DIR_TMP_CDF ; fi
if [ ! -d $DIR_TMP_LOG  ] ; then mkdir $DIR_TMP_LOG ; fi
if [ ! -d $DIR_TMP_RUN  ] ; then mkdir $DIR_TMP_RUN ; fi

## ======================================================================
## FILES NAMES
## ======================================================================
EXEC_NAME=nemo.exe
BATMET_FIL=bathy_salishsea_1500m_20210406.nc
COORD_FIL=coordinates_salishsea_1500m.nc
BATLE_FIL=
BFRI_FIL=
SHLAT_FIL=
# INITIAL CONDITIONS FILES 
TEMP_FIL=
SAL_FIL=
U_FIL=
V_FIL=
SSH_FIL=
# RIVER RUNOFF FILES                                            #jpp - names must be consistent with
RUNOFF_FIL=
RUNOFF_TEMP_FIL=
RUNOFF_MASK_FIL=
# CHLOROPHYL related files
CHLA_FIL=
RGB_FIL=
# SEA ICE FILES
INITICE_FIL=''							#jpp- not used
CICE_GRID='' 							#jpp- not used
# TIDAL FILES
compo_list='M2'                                                 #!!!JPP TO BE CLEANED !!!
TIDE_PREF=							#jpp- not used, see comment top
TIDEH_FIL=
TIDEU_FIL=
TIDEV_FIL=
TMX_K1=                                                         #jpp- not used
TMX_M2=                                                         #jpp- not used
TMX_MASK=                                                       #jpp- not used
TMX_MASK=     
# PISCES FILES 
DUST_FIL=
RIVER_FIL=
# OBC/BDY files 
obc_list=
obcpref2D=
obcsuf2D=
obcprefNorth=
                                                                #jpp- leave empty (i.e. '') not to look for files
# Weights files for atmospheric forcing
WEIGHTS_FIL=weights_bilin_OPPwest.nc

RESTART_FILE=${HOME}/MEOPAR/forcing/restart/hindcast.201812/31dec15/SalishSea_01682640_restart.nc

## ======================================================================
## FUNCTIONS rappatrie, rappatrest and expatrie
## ======================================================================
	rappatrie ()
	{
	    if [ -f $2/$1 ] ; then
		ln -sf $2/$1 $1 ;
	    else
		cp $3 $1;
		cp $1 $2/$1 ;
	    fi ;
	}
	rappatrie2 ()
	{
	    if [ -f $3/$2 ] ; then
		echo  $3/$2 PRESENT;
		ln -sf $3/$2 $1 ;
	    else
		echo  $3/$2 NOT PRESENT ;
		tarfile=ERAinterim_$6_$5.tar ;
		cp $4/${tarfile}.gz $3 ;
		cd $3 ;
		gunzip $3/${tarfile}.gz ;
		tar xvf $3/${tarfile} ;
		cd ${DIR_TMP_RUN} ; 
		ln -sf  $3/$2 $1 ;
	    fi ;
	}
	rappatrest ()
	{
	    if [ -f $2 ] ; then
		ln -sf $2 $1 ;
	    else
		cp $3 $1 ;
	    fi ;
	}
        expatrie () {
            cp $2/$1 $3/$1 ;
	}
	expatrie_mv () {		 # JPP - move faster than copy
	    if [ ! -d $3 ] ; then
		echo 'mkdir directory: '$3
		mkdir -m 755 -p $3
	    fi
	    mv $2/$1 $3/$1 ;
	}

## ======================================================================
## FUNCTION j2d
## ======================================================================
	j2d ()
	{
#            set -axe
#	    typeset -i day month year tmpday centuries
#	    ((tmpday = $1 + 712164))
#	    ((centuries = (4 * tmpday - 1) / 146097))
#	    ((tmpday += (centuries - centuries/4)))
#	    ((year = (4 * tmpday - 1) / 1461))
#	    ((tmpday -= (1461 * year) / 4))
#	    ((month = (10 * tmpday - 5) / 306))
#	    ((day = tmpday - (306 * month + 5) / 10))
#	    ((month += 2))
#	    ((year += month/12))
#	    ((month = month % 12 + 1))
            set -axe
#           typeset -i day month year tmpday centuries
           let tmpday=$(($1+712164))
           let centuries=$(((4*$tmpday-1)/146097))
           let tmpday=$(($tmpday + ($centuries - $centuries/4)))
           let year=$(( (4 * $tmpday - 1) / 1461))
	   let tmpday=$(( $tmpday- (1461 * $year) / 4))
           let month=$(( (10 * $tmpday - 5) / 306))
           let day=$(( $tmpday - (306 * $month + 5) / 10))
           let month=$(( $month + 2))
           let year=$(( $year + $month/12))
           let month=$(( $month % 12 + 1))
	    echo ${year}`printf %2.2i ${month}``printf %2.2i ${day}`
	}

## ======================================================================
## FUNCTION d2j
## ======================================================================
	d2j ()
	{
	    typeset -i day month year tmpmonth tmpyear
	    ndate=$1
	    year=`echo ${ndate} | cut -c 1-4`
	    month=`echo ${ndate} | cut -c 5-6`
	    day=`echo ${ndate} | cut -c 7-8`
	    ((tmpmonth = 12 * year + month - 3))
	    ((tmpyear = tmpmonth / 12))
	    echo $(((734*tmpmonth + 15)/24 - 2*tmpyear + tmpyear/4 - tmpyear/100 + tmpyear/400 + day - 712164))
	}

## ======================================================================
## Print informations relative to the run
## ======================================================================
#. ./$BACKEND

set +x
echo
echo Configuration name : $CONFIG
echo Number of procs : $NBPROC
echo
if [ $RESTART == "yes" ] ; then
    echo Hot restart from date $NDATEDEB
else
    echo Initialisation from climatology or analytical profiles
fi
echo
echo
echo home : $HOME
echo 
echo "CALCU directories :"
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo $DIR_CALCU_EXE    : directory where is the NEMO executable
echo $DIR_CALCU_SHR    : directory where is the shared files
echo $DIR_CALCU_XIO    : directory where is the XIOS server
echo $DIR_CALCU_CTL    : directory where are control files of the run
echo $DIR_CALCU_UTL    : directory where to find utilities
echo
echo "MACHINEARCHI" directories :
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo $DIR_ARCHI_ATMOS      : directory to stock the atmos forcing fields
echo $DIR_ARCHI_INPUT     : directory for the input files
echo $DIR_ARCHI_CLIM      : directory for the climatology files
echo $DIR_ARCHI_OBC       : directory for the obc files
echo $DIR_ARCHI_HOME      : directory to stock the results
echo $DIR_ARCHI_RST       : directory to stock the restarts
echo $DIR_ARCHI_CDF       : directory for combined netCDF results
echo $DIR_ARCHI_LOG       : directory for ice something ...
echo
echo "TEMP" directories :
echo ~~~~~~~~~~~~~~~~~~~~
echo $DIR_TMP_HOME   : home temp directory 
echo $DIR_TMP_RST    : directory for the restart
echo $DIR_TMP_CDF    : directory for the cdf
echo $DIR_TMP_LOG    : directory for the log
echo $DIR_TMP_RUN    : directory to run
echo $DIR_TMP_ATM  : directory for the atmos forcing
echo $DIR_TMP_OBC    : directory for the obc files
echo
echo Run features :
echo ~~~~~~~~~~~~~~
echo Executable name : $EXEC_NAME
echo Bathymetry : $BATMET_FIL
echo Coordinates : $COORD_FIL
echo K1 tidal mixing : $TMX_K1
echo M2 tidal mixing : $TMX_M2
echo Mask tidal mixing : $TMX_MASK
echo Shlat : $SHLAT_FIL
echo Clim T : $TEMP_FIL
echo Clim S : $SAL_FIL
echo Clim U : $U_FIL
echo Clim V : $V_FIL
echo Runoff : $RUNOFF_FIL
echo Chlorophyl : $CHLA_FIL
echo List of the open boundaries : $obc_list
echo
