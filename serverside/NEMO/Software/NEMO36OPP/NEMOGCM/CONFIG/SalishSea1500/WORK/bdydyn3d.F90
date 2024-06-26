MODULE bdydyn3d
   !!======================================================================
   !!                       ***  MODULE  bdydyn3d  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on baroclinic velocities
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite 
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
#if defined key_opp
   !!            3.6  !  2017     (S. Allen) sponge layer
#endif
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   !!   bdy_dyn3d        : apply open boundary conditions to baroclinic velocities
   !!   bdy_dyn3d_frs    : apply Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdylib          ! for orlanski library routines
#if defined key_opp
   USE ldfdyn_oce      ! for sponge layer
#endif
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  !
   Use phycst

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn3d     ! routine called by bdy_dyn
   PUBLIC   bdy_dyn3d_dmp ! routine called by step

   !! * Substitutions
#  include "domzgr_substitute.h90"
#if defined key_opp
#  include "ldfdyn_substitute.h90"
#endif
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdydyn3d.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn3d( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for baroclinic velocities
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt     ! Main time step counter
      !!
      INTEGER               :: ib_bdy ! loop index
      !!
#if defined key_opp
      IF( kt == nit000 ) THEN
        DO ib_bdy=1, nb_bdy
          IF (cn_dyn3d(ib_bdy) == 'orlanski_w_sponge') THEN
            CALL bdy_dyn3d_sponge( idx_bdy(ib_bdy), ib_bdy )
          ELSEIF (cn_dyn3d(ib_bdy) == 'sponge') THEN
            CALL bdy_dyn3d_sponge( idx_bdy(ib_bdy), ib_bdy )
            cn_dyn3d(ib_bdy) = 'none'
          ENDIF
        ENDDO
      ENDIF
#endif
      DO ib_bdy=1, nb_bdy

         SELECT CASE( cn_dyn3d(ib_bdy) )
         CASE('none')
            CYCLE
         CASE('frs')
            CALL bdy_dyn3d_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('specified')
            CALL bdy_dyn3d_spe( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('zero')
            CALL bdy_dyn3d_zro( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
#if defined key_opp
         CASE('orlanski')
            CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.false., ll_spg=.false.)
         CASE('orlanski_npo')
            CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.true. , ll_spg=.false.)
         CASE('orlanski_w_sponge')
            CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.false., ll_spg=.true. )
#else
         CASE('orlanski')
            CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.false. )
         CASE('orlanski_npo')
            CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.true. )
#endif
         CASE DEFAULT
            CALL ctl_stop( 'bdy_dyn3d : unrecognised option for open boundaries for baroclinic velocities' )
         END SELECT
      ENDDO

   END SUBROUTINE bdy_dyn3d

   SUBROUTINE bdy_dyn3d_spe( idx, dta, kt , ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_spe  ***
      !!
      !! ** Purpose : - Apply a specified value for baroclinic velocities
      !!                at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_spe')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            ua(ii,ij,jk) = dta%u3d(jb,jk) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            va(ii,ij,jk) = dta%v3d(jb,jk) * vmask(ii,ij,jk)
         END DO
      END DO
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ! Boundary points should be updated  
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_spe')

   END SUBROUTINE bdy_dyn3d_spe

   SUBROUTINE bdy_dyn3d_zro( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_zro  ***
      !!
      !! ** Purpose : - baroclinic velocities = 0. at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   ib, ik         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd, zcoef   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_zro')
      !
      igrd = 2                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ua(ii,ij,ik) = 0._wp
         END DO
      END DO

      igrd = 3                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            va(ii,ij,ik) = 0._wp
         END DO
      END DO
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ;   CALL lbc_bdy_lnk( va, 'V', -1.,ib_bdy )   ! Boundary points should be updated
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_zro')

   END SUBROUTINE bdy_dyn3d_zro

   SUBROUTINE bdy_dyn3d_frs( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for baroclinic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_frs')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta%u3d(jb,jk) - ua(ii,ij,jk) ) ) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta%v3d(jb,jk) - va(ii,ij,jk) ) ) * vmask(ii,ij,jk)
         END DO
      END DO 
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_frs')

   END SUBROUTINE bdy_dyn3d_frs


#if defined key_opp
   SUBROUTINE bdy_dyn3d_orlanski( idx, dta, ib_bdy, ll_npo, ll_spg )
#else
   SUBROUTINE bdy_dyn3d_orlanski( idx, dta, ib_bdy, ll_npo )
#endif
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn3d_orlanski  ***
      !!             
      !!              - Apply Orlanski radiation to baroclinic velocities. 
      !!              - Wrapper routine for bdy_orlanski_3d
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! BDY set index
      LOGICAL,                      INTENT(in) ::   ll_npo  ! switch for NPO version
#if defined key_opp
      LOGICAL,                      INTENT(in) ::   ll_spg  ! switch for SPONGE version
#endif

      INTEGER  ::   jb, igrd                               ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_orlanski')
      !
      !! Note that at this stage the ub and ua arrays contain the baroclinic velocities. 
      !
      igrd = 2      ! Orlanski bc on u-velocity; 
      !            
#if defined key_opp
      CALL bdy_orlanski_3d( idx, igrd, ub, ua, dta%u3d, ll_npo, ll_spg )
#else
      CALL bdy_orlanski_3d( idx, igrd, ub, ua, dta%u3d, ll_npo )
#endif

      igrd = 3      ! Orlanski bc on v-velocity
      !  
#if defined key_opp
      CALL bdy_orlanski_3d( idx, igrd, vb, va, dta%v3d, ll_npo, ll_spg )
#else
      CALL bdy_orlanski_3d( idx, igrd, vb, va, dta%v3d, ll_npo )
#endif
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_orlanski')
      !
   END SUBROUTINE bdy_dyn3d_orlanski


   SUBROUTINE bdy_dyn3d_dmp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_dmp  ***
      !!
      !! ** Purpose : Apply damping for baroclinic velocities at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::  ib_bdy          ! loop index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_dmp')
      !
      !-------------------------------------------------------

      DO ib_bdy=1, nb_bdy
         IF ( ln_dyn3d_dmp(ib_bdy) .and. cn_dyn3d(ib_bdy) /= 'none' ) THEN
            igrd = 2                      ! Relaxation of zonal velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%u3d(jb,jk) - &
                                   ub(ii,ij,jk) + ub_b(ii,ij)) ) * umask(ii,ij,jk)
               END DO
            END DO
            !
            igrd = 3                      ! Relaxation of meridional velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%v3d(jb,jk) -  &
                                   vb(ii,ij,jk) + vb_b(ii,ij)) ) * vmask(ii,ij,jk)
               END DO
            END DO
         ENDIF
      ENDDO
      !
      CALL lbc_lnk( ua, 'U', -1. )   ;   CALL lbc_lnk( va, 'V', -1. )   ! Boundary points should be updated
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_dmp')

   END SUBROUTINE bdy_dyn3d_dmp


#if defined key_opp
   SUBROUTINE bdy_dyn3d_sponge( idx, ib_bdy)
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_sponge  ***
      !!
      !! ** Purpose : Increase horizontal viscosity at boundaries
      !!
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      INTEGER, POINTER :: nbr      ! short cut
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_sponge')
      !
      IF( ln_dynldf_iso) CALL ctl_stop &
      & ( 'bdy_dyn3d : Isoneutral viscosity uses only a single diffusivity value, cannot use sponge' )
      igrd = 1                     ! viscosity defined on T-grid
      DO jb = 1, idx%nblen(igrd)
#if defined key_dynldf_c3d
            DO jk = 1, jpkm1
               ii   = idx%nbi(jb,igrd)
               ij   = idx%nbj(jb,igrd)
	       nbr => idx%nbr(jb,igrd)
	       zwgt = rn_max_sponge(ib_bdy)*cos(rpi*(nbr-1.)/(2.*(nn_rimwidth(ib_bdy))))
	       fsahmt(ii, ij, jk) = fsahmt(ii, ij, jk) + zwgt
            END DO
#elif defined key_dynldf_c2d
	 jk = 1
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
	 nbr => idx%nbr(jb,igrd)
	 zwgt = rn_max_sponge(ib_bdy)*cos(rpi*(nbr-1.)/(2.*(nn_rimwidth(ib_bdy))))
	 fsahmt(ii, ij, jk) = fsahmt(ii, ij, jk) + zwgt
#endif
      END DO
      DO jb = 1, idx%nblen(igrd)
#if defined key_dynldf_c3d
	    DO jk = 1, jpkm1
               ii   = idx%nbi(jb,igrd)
               ij   = idx%nbj(jb,igrd)
	       fsahmf(ii, ij, jk) = 0.25 * (fsahmt(ii, ij, jk) + fsahmt(ii+1, ij, jk) + &
                                         fsahmt(ii, ij+1, jk) + fsahmt(ii+1, ij+1, jk))
            END DO
#elif defined key_dynldf_c2d
	 jk = 1
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
	 fsahmf(ii, ij, jk) = 0.25 * (fsahmt(ii, ij, jk) + fsahmt(ii+1, ij, jk) + &
                                         fsahmt(ii, ij+1, jk) + fsahmt(ii+1, ij+1, jk))
#endif
      END DO

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_sponge')

   END SUBROUTINE bdy_dyn3d_sponge
#endif


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn3d( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d

   SUBROUTINE bdy_dyn3d_dmp( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d_dmp

#endif

   !!======================================================================
END MODULE bdydyn3d
