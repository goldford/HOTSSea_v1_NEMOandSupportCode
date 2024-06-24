#!/bin/ksh

set -ex

## ======================================================================
## AUTO-RESTART: POST-PROCESSING SCRIPT
## POST-TREATMENT.ksh
## ======================================================================
# Change to control directory if given as command line parameter by ord_soumet
if [ ! -z "$1" ] && [ -d $1 ]; then cd $1; fi

# Import run settings
. ./environment.ksh

## ======================================================================
## Read run number, nit000 and nitend from last line of $CONFIG.db
## ======================================================================
read no nit000 nitend < <(tail -1 $CONFIG.db)
typeset -Z8 nit000 nitend

## ======================================================================
## Rebuild restart files
## ======================================================================
cd $DIR_TMP_RST
for restfile in ${CONFIG}_*_restart_0000.nc; do
  rebuild_successful=0
  if [ ! -x ${DIR_CALCU_RBL}/BLD/bin/rebuild_nemo.exe ] ; then break; fi
  if [ ! -f ${restfile} ]; then continue; fi
  nitrst=$restfile
  typeset -R24 nitrst  # trim rightmost 24 characters
  typeset -L8  nitrst  # trim leftmost 8 characters
  typeset -Z8  nitrst  # ensure 8 digits with leading zeros
  nombase=${CONFIG}_${nitrst}_restart
  echo "Rebuild lddcheck"
  hostname
  echo ldd ${DIR_CALCU_RBL}/BLD/bin/rebuild_nemo.exe
  ldd ${DIR_CALCU_RBL}/BLD/bin/rebuild_nemo.exe
  echo "Rebuilding $nombase"
  $DIR_CALCU_RBL/rebuild_nemo -t 1 $nombase $NBPROC > rebuild.log 2>&1
  rebuild_successful=`grep "NEMO rebuild completed successfully" rebuild.log | wc -l`
  if [ $rebuild_successful -eq 1 ]; then
    echo "Restart rebuild of $nombase successful"
    # cleanup per-proc files
    rm -f rebuild.log nam_rebuild
    typeset -Z4 jproc
    for jproc in `seq 0 $(($NBPROC-1))`; do
      nomrestart=${CONFIG}_${nitrst}_restart_${jproc}.nc
      rm -f $nomrestart
    done
  fi
  unset nitrst
done

## ======================================================================
## Rebuild output.init and mesh_mask if exists
## ======================================================================
cd $DIR_TMP_CDF
for prefix in mesh_mask output.init output.abort; do
  if [ ! -x ${DIR_CALCU_RBL}/BLD/bin/rebuild_nemo.exe ] ; then break; fi
  if [ ! -f ${DIR_TMP_CDF}/${prefix}_0000.nc ]; then continue; fi
  rebuild_successful=0
  $DIR_CALCU_RBL/rebuild_nemo -t 1 $prefix $NBPROC > rebuild.log
  rebuild_successful=`grep "NEMO rebuild completed successfully" rebuild.log | wc -l`
  if [ $rebuild_successful -eq 1 ]; then
    nom=${prefix}.nc
    nomdest=${prefix}_${nit000}.nc
    echo "Rebuild for $nom successful, renaming it to $nomdest"
    mv $nom $nomdest
    # cleanup per-proc files
    rm -f rebuild.log nam_rebuild
    typeset -Z4 jproc
    for jproc in `seq 0 $(($NBPROC-1))`; do
      rm -f ${prefix}_${jproc}.nc
    done
  fi
done

## ======================================================================
## Special case renaming salish sea files
## ======================================================================
set +e
#if [ "$CONF_NAME" == "SalishSea" ] ; then
#  for fname in [A-Z]*.nc; do
#    if [[ $fname != ${CONF_NAME}* ]]; then
#      nombase=`basename $fname .nc`
#      mv $fname ${nombase}_${nit000}_${nitend}.nc
#    fi
#  done
#fi

## ======================================================================
## Move results to archive directories
## ======================================================================
rsync -axvH --no-g --no-p ${DIR_TMP_LOG}/* $DIR_ARCHI_LOG && rm ${DIR_TMP_LOG}/*
rsync -axvH --no-g --no-p ${DIR_TMP_CDF}/* $DIR_ARCHI_CDF && rm ${DIR_TMP_CDF}/*
rsync -axvH --no-g --no-p ${DIR_TMP_RST}/* $DIR_ARCHI_RST && rm ${DIR_TMP_RST}/*
#mv -f ${DIR_TMP_LOG}/* $DIR_ARCHI_LOG
#mv -f ${DIR_TMP_CDF}/* $DIR_ARCHI_CDF
#mv -f ${DIR_TMP_RST}/* $DIR_ARCHI_RST
set -e

## ======================================================================
## Automatic re-submit by calling RAPPAT.ksh
## ======================================================================
cd $DIR_CALCU_CTL
if [ $no -lt $MAXSUB ] ; then
  # Add a new line to the database and resubmit via RAPPAT.ksh
  typeset +Z nit000 nitend
  let no=$(($no+1))
  let dif=$(($nitend-$nit000+1))
  let nit000=$(($nitend+1))
  let nitend=$(($nitend+dif))
  echo $no $nit000 $nitend >> $CONFIG.db
  echo "Automatic resubmit for job $no"
  ksh ./RAPPAT.ksh
fi

