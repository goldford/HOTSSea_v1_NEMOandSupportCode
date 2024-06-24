#!/bin/ksh

set -ex

## ======================================================================
## AUTO-RESTART: NEMO RUN SCRIPT
## $CONF.ksh
## ======================================================================
# Change to control directory given as command line parameter by ord_soumet
if [ ! -z "$1" ] && [ -d $1 ]; then cd $1; fi

# Import run settings
. ./environment.ksh

## ======================================================================
## XIOS Server (deteched mode) 
## ======================================================================
if [ $NBPROC_XIO -gt 0 ] ; then
    XIOS=': -np '${NBPROC_XIO}' ./xios_server.exe'
else
    XIOS=''
fi

## ======================================================================
## Execution of OPA
## ======================================================================
echo '*******************'
echo '* RUN * RUN * RUN *'
echo '*******************'

cd $DIR_TMP_RUN
pwd

time mpirun -np ${NBPROC} ./opa ${XIOS}  | tee > nemo.stdout 2> nemo.stderr

## ======================================================================
## Read run number, nit000 and nitend from last line of $CONFIG.db
## ======================================================================
read no nit000 nitend < <(tail -1 $DIR_CALCU_CTL/$CONFIG.db)

## ======================================================================
## Check if run was successful
## ======================================================================
run_success=0
if [ -f time.step ]; then
  read nitstp < time.step
  if [ $nitstp -eq $nitend ] ; then
    run_success=1
  fi
fi

## ======================================================================
## Move files we wish to keep from RUN dir to LOG/OUT/RST dirs
## ======================================================================
set +e
typeset -Z3 no
typeset -Z8 nit000 nitend
for fname in *timing.output *ocean.output *time.step *layout.dat *namelist_??? *solver.stat *output.namelist.dyn nemo.stdout nemo.stderr; do
  mv ${fname} ${DIR_TMP_LOG}/${fname}_R${no}_${nit000}_${nitend}
done

## Move Restarts and outputs out of the execution directory
## (Archives and definitive move done in POST_TREATMENT.ksh)
if [ "$CONF_NAME" == "SalishSea" ] ; then
  mv [A-Z]*.nc        ${DIR_TMP_CDF}  #- SalishSea specifc outputs
fi
mv ${CONFIG}*.nc      ${DIR_TMP_CDF}  #- Model output
mv output.init*.nc    ${DIR_TMP_CDF}  #- Initial conditions if active
mv output.abort*.nc   ${DIR_TMP_CDF}  #- Abort file
mv mesh_mask*.nc      ${DIR_TMP_CDF}  #- Mesh mask if active
mv OUTRST/*restart*nc ${DIR_TMP_RST}  #- Restart files
set -e

## ======================================================================
## Submit job for post processing with POST-TREATMENT.ksh
## ======================================================================

if [ $run_success -eq 0 ]; then
  exit 1
fi

#echo "Submitting POST-TREATMENT.ksh job"
#cd $DIR_CALCU_CTL
#ord_soumet ./POST-TREATMENT.ksh -image ${CONTAINER} -mach ${MACHINE} -cpus 1 -mpi -cm ${MEMLIMIT} \
 #                               -t ${TIMELIMIT_PP} -listing ${DIR_CALCU_CTL} -q ${QUEUE} \
 #                               -o ${DIR_CALCU_CTL} -project ${PROJECT} -jn ${CONF_NBR}-PP \
 #                               -mail Maxim.Krassovski@dfo-mpo.gc.ca -notify -keep

## ======================================================================
## End of the $CONF script
## ======================================================================
