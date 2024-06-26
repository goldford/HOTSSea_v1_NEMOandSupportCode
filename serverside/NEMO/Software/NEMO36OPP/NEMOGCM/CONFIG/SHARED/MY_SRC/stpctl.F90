MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! surface boundary conditions variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE lib_fortran     ! Fortran routines library 
   USE dynspg_oce      ! pressure gradient schemes 
   USE c1d             ! 1D vertical configuration


   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: stpctl.F90 7847 2017-03-30 13:30:29Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol 
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !! ** Actions :   'time.step' file containing the last ocean time-step
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence
      !!
      CHARACTER(len = 32) ::        clfname ! time stepping output file name
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ii, ij, ik              ! temporary integers
      REAL(wp) ::   zumax, zsmin, zssh2, zsshmax     ! temporary scalars
      INTEGER, DIMENSION(3) ::   ilocu      ! 
      INTEGER, DIMENSION(2) ::   ilocs      ! 
      !!----------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         ! open time.step file with special treatment for SAS
         IF ( nn_components == jp_iam_sas ) THEN
            clfname = 'time.step.sas'
         ELSE
            clfname = 'time.step'
         ENDIF
         CALL ctl_opn( numstp, TRIM(clfname), 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
      ENDIF

      IF(lwp) WRITE ( numstp, '(1x, i8)' )   kt      !* save the current time step in numstp
      IF(lwp) REWIND( numstp )                       !  --------------------------

#if defined key_opp
      ! Abort simulation if NaN found in SSH (code courtesy Susan Allen)
      !                                              ! * Test SSH !
      zssh2 = SUM( sshn(:,:) * sshn(:,:) * tmask_i(:,:) )
      IF( lk_mpp )   CALL mpp_sum( zssh2 )      ! sum over the global domain

      IF (isnan(zssh2)) THEN
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl: sea surface height is nan'
            WRITE(numout,*) ' ====== '
            WRITE(numout,9401) kt
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         kindic = -3
      ENDIF
9401  FORMAT (' kt=',i6,' Nan found in SSH.')
#endif
      !                                              !* Test maximum of velocity (zonal only)
      !                                              !  ------------------------
      !! zumax = MAXVAL( ABS( un(:,:,:) ) )                ! slower than the following loop on NEC SX5
      zumax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zumax = MAX(zumax,ABS(un(ji,jj,jk)))
          END DO 
        END DO 
      END DO        
      IF( lk_mpp )   CALL mpp_max( zumax )                 ! max over the global domain
      !
      IF( MOD( kt, nwrite ) == 1 .AND. lwp )   WRITE(numout,*) ' ==>> time-step= ',kt,' abs(U) max: ', zumax
      !
      IF( zumax > 20.e0 ) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(un),umask,zumax,ii,ij,ik)
         ELSE
            ilocu = MAXLOC( ABS( un(:,:,:) ) )
            ii = ilocu(1) + nimpp - 1
            ij = ilocu(2) + njmpp - 1
            ik = ilocu(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl: the zonal velocity is larger than 20 m/s'
            WRITE(numout,*) ' ====== '
            WRITE(numout,9400) kt, zumax, ii, ij, ik
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         kindic = -3
      ENDIF
9400  FORMAT (' kt=',i6,' max abs(U): ',1pg11.4,', i j k: ',3i5)

      !                                              !* Test minimum of salinity
      !                                              !  ------------------------
      !! zsmin = MINVAL( tsn(:,:,1,jp_sal), mask = tmask(:,:,1) == 1.e0 )  slower than the following loop on NEC SX5
      zsmin = 100.e0
      DO jj = 2, jpjm1
         DO ji = 1, jpi
            IF( tmask(ji,jj,1) == 1) zsmin = MIN(zsmin,tsn(ji,jj,1,jp_sal))
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_min( zsmin )                ! min over the global domain
      !
      IF( MOD( kt, nwrite ) == 1 .AND. lwp )   WRITE(numout,*) ' ==>> time-step= ',kt,' SSS min:', zsmin
      !
      IF( zsmin < 0.) THEN 
         IF (lk_mpp) THEN
            CALL mpp_minloc ( tsn(:,:,1,jp_sal),tmask(:,:,1), zsmin, ii,ij )
         ELSE
            ilocs = MINLOC( tsn(:,:,1,jp_sal), mask = tmask(:,:,1) == 1.e0 )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         ENDIF
         !
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) 'stp_ctl : NEGATIVE sea surface salinity'
            WRITE(numout,*) '======= '
            WRITE(numout,9500) kt, zsmin, ii, ij
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         kindic = -3
      ENDIF
9500  FORMAT (' kt=',i6,' min SSS: ',1pg11.4,', i j: ',2i5)

      
      IF( lk_c1d )  RETURN          ! No log file in case of 1D vertical configuration

      ! log file (solver or ssh statistics)
      ! --------
      IF( lk_dynspg_flt ) THEN      ! elliptic solver statistics (if required)
         !
         IF(lwp) WRITE(numsol,9200) kt, niter, res, SQRT(epsr)/eps       ! Solver
         !
         IF( kindic < 0 .AND. zsmin > 0.e0 .AND. zumax <= 20.e0 ) THEN   ! create a abort file if problem found 
            IF(lwp) THEN
               WRITE(numout,*) ' stpctl: the elliptic solver DO not converge or explode'
               WRITE(numout,*) ' ====== '
               WRITE(numout,9200) kt, niter, res, sqrt(epsr)/eps
               WRITE(numout,*)
               WRITE(numout,*) ' stpctl: output of last fields'
               WRITE(numout,*) ' ======  '
            ENDIF
         ENDIF
         !
      ELSE                                   !* ssh statistics (and others...)
         IF( kt == nit000 .AND. lwp ) THEN   ! open ssh statistics file (put in solver.stat file)
            CALL ctl_opn( numsol, 'solver.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         ENDIF
         !
         zsshmax = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,1) == 1) zsshmax = MAX( zsshmax, ABS(sshn(ji,jj)) )
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_max( zsshmax )                ! min over the global domain
         !
         IF( MOD( kt, nwrite ) == 1 .AND. lwp )   WRITE(numout,*) ' ==>> time-step= ',kt,' ssh max:', zsshmax
           !
           IF( zsshmax > 15.e0 ) THEN
              IF (lk_mpp) THEN
                CALL mpp_maxloc( ABS(sshn(:,:)),tmask(:,:,1),zsshmax,ii,ij)
              ELSE
                ilocs = MAXLOC( ABS(sshn(:,:)) )
                ii = ilocs(1) + nimpp - 1
                ij = ilocs(2) + njmpp - 1
              ENDIF
              !
              IF(lwp) THEN
                 WRITE(numout,cform_err)
                 WRITE(numout,*) 'stp_ctl : the ssh is larger than 15m'
                 WRITE(numout,*) '======= '
                 WRITE(numout,9600) kt, zsshmax, ii, ij
                 WRITE(numout,*)
                 WRITE(numout,*) '          output of last fields in numwso'
              ENDIF
             kindic = -3
         ENDIF
9600     FORMAT (' kt=',i6,' max ssh: ',1pg11.4,', i j: ',2i5)

         zssh2 = glob_sum( sshn(:,:) * sshn(:,:) )
         !
         IF(lwp) WRITE(numsol,9300) kt, zssh2, zumax, zsmin      ! ssh statistics
         !
      ENDIF
9200  FORMAT('it:', i8, ' iter:', i4, ' r: ',d23.16, ' b: ',d23.16 )
9300  FORMAT(' it :', i8, ' ssh2: ', d23.16, ' Umax: ',d23.16,' Smin: ',d23.16)
      !
   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
