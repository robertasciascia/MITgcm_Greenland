C $Header: /u/gcmpack/MITgcm_contrib/icefront/3D_example/code/icefront_thermodynamics.F,v 1.1 2014/01/24 23:03:22 RS $
C $Name:  $

CBOP
C     !SMR_PARAMETERIZATION
C     !INTERFACE:
    SUBROUTINE SMR_PARAMETERIZATION(
    I                        myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *=============================================================*
C     | S/R  SMR_PARAMETERIZATION
C     |
C     |   Routine that given a runoff value computes sumbarine
C     |   melting and the density of the glacially modified waters
C     |   at neutral buoyancy (after plume dynamics)
C     |
C     *=============================================================*

C     !USES:
    IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myIter :: iteration counter for this thread
C     myTime :: time counter for this thread
C     myThid :: thread number for this instance of the routine.
    _RL  myTime
    INTEGER myIter
    INTEGER myThid

CEOP

C     !LOCAL VARIABLES :
C     === Local variables ===
C     I,J,K,bi,bj,pk        :: loop counters
C     tLoc, rLoc            :: local in-situ density and sum of in situ temperature
C     smr                   :: sumbarine melting
C     qsg                   :: subglacial discharge
C     tMean                 :: mean temperature of the water column
C     rhoP, tempP, saltP    :: plume temperature, salinity and density @ z
C     rhoPo, tempPo, saltPo :: plume temperature, salinity and density @ z-1
C     QPo                   :: plume flow rate @ z-1
C     g_prime               :: reduced gravity 
C     B, gP,wP,QP           :: Plume properties @ z
C     a0, b0, c0, d0        :: Coefficients to compute smr from qsg
C     tempA, saltA          :: In-situ tempererature and salinity


      INTEGER I,J,K,pk
      INTEGER bi,bj
      _RL qsg_fj(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL rho_fj(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL temp_fj(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL salt_fj(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL rhoP, tempP, saltP
      _RL QP
      _RL tLoc, rLoc
      _RL tempA, saltA
      _RL g_prime
      _RL a0, b0, c0,d0
      _RL tMean
      _RL qsg, smr
      _RL g_prime, B, gP,wP
      _RL tempPo, saltPo, rhoPo, Qo


     _RL SW_TEMP
     EXTERNAL SW_TEMP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


# ifdef ALLOW_ADDFLUID
C     Parameters to compute smr from runoff qsg in m^3/s

      a0 = 0.01      _d  0
      b0 = 0.65      _d  0
      c0 = 0.002     _d  0
      d0 = 1.35      _d  0
C     Initialize variables
      B = 0.         _d  0
      g_prime = 0.   _d  0
      gP = 0.        _d  0
      wP = 0.        _d  0
      rhoP = 0.      _d  0
      QP = 0.        _d  0
      tempP = 0.     _d  0
      saltP = 0.     _d  0
      rhoPo = 0.     _d  0
      QPo = 0.       _d  0
      tempPo = 0.    _d  0
      saltPo = 0.    _d  0
      qsg_fj= 0.     _d  0
      rho_fj=0.      _d  0
      temp_fj=0.     _d  0
      salt_fj=0.     _d  0

C     Do the calculation only if Addmass GT 0

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

         DO J = 1, sNy
          DO I = 1, sNx
             k=kLowC(i,j,bi,bj)

             IF(addMass(I,J,k,bi,bj) .GT. 0. _d 0) THEN
                tLoc=sum(theta(I,J,1:Nr,bi,bj),3)                 ! compute the depth averaged temperaure
                tMean=tLoc/k                                      !
                qsg=addMass(I,J,k,bi,bj)/rhoConstFresh            ! Addmass in kg/s convert to m^3/s
                smr= ((a0*(qsg)**b0)+c0)*(tMean**d0)              ! compute smr
                qsg=qsg+smr                                       ! add SMR to QSG

C From Bush JFM 99 and Turner book LINE plume models to compute the intrusion level
C z=0
                qsg=qsg/delY(J);                                  ! flow rate per unit length
                rLoc=rhoanoma(I,J,k,bi,bj)+rhoConst               ! ambient water density
                g_prime=gravity*(rLoc-rhoConstFresh)/rLoc         ! reduced gravity
                B=2*g_prime*qsg;                                  ! buoyancy flux per unit length and
                                                                  ! multiplied by two (against the glacier)

                gP=(2.6*(B**(2/3)))./(k*delZ(k));                 ! plume buoyancy @ z
                wP=1.8*(B**(1/3));                                ! plume velocity @ z
                QP(=(B/gP);                                       ! plume flow rate in m^2/s @ z
                tempA=theta(I,J,k,bi,bj);                         ! ambient temperature @ z
                saltA=salt(I,J,k,bi,bj);                          ! ambiente salinity @ z
                rhoP=((rhoConstFresh*qsg)+(rLoc*(QP-qsg)))/QP     ! plume density @ z
                tempP=((temp_addMass*qsg)+(tempA*(QP-qsg)))/QP    ! plume temperature @ z
                saltP=((salt_addMass*qsg)+(saltA*(QP-qsg)))/QP    ! plume salinity @ z

C start loop over z (from seafloor to the free surface )

                DO pk = k-1,1,-1                                  ! compute plume bouancy and density

                    IF((rLoc-rhoP)) .GT. 0.01) THEN               ! compare plume density with ambient density

                        rhoPo=rhoP                                ! store previous densitity @ z-1
                        tempPo=tempP                              ! store previous temperature @ z-1
                        saltPo=saltP                              ! store previous velocity @ z-1
                        QPo=QP                                    ! store previous flow fate @ z-1
                        rLoc=rhoanoma(I,J,pk,bi,bj)+rhoConst      ! update local density @ z
                        gP=(2.6*(B**(2/3)))./(pk*delZ(pk));       ! plume buoyancy @ z
                        wP=1.8*(B**(1/3));                        ! plume velocity @ z

                        QP(=(B/gP);                               ! plume flow rate in m^2/s @ z
                        tempA=theta(I,J,pk,bi,bj);                ! ambient temperature @ z
                        satlA=salt(I,J,pk,bi,bj);                 ! ambient salinity    @ z
                        rhoP=((rhoPo*QPo)+(rLoc*(QP-QPo)))/QP     ! plume density       @ z
                        tempP=((tempPo*QPo)+(rLoc*(QP-QPo)))/QP   ! plume temperature   @ z
                        saltP=((saltPo*QPo)+(rLoc*(QP-QPo)))/QP   ! plume salinity      @ z


                    ENDIF                                         ! if (rLoc-rhoP)) .GT. 0.01
                ENDDO                                             ! pk = k,1

C final quantities to be passed again to the addmass field  
                qsg_fj(I,J,pk,bi,bj)=(QP/2)*rhoConstFresh         ! modified qsg (qsg+smr+entrainment)

                rho_fj(I,J,pk,bi,bj)=rhoP                         ! modified qsg density 
                temp_fj(I,J,pk,bi,bj)=tempP                       ! modified qsg temperature
                salt_fj(I,J,pk,bi,bj)=saltP                       ! modified qsg salinity

             ENDIF              ! addmass GT 0


          ENDDO                 ! I = 1, sNx
         ENDDO                  ! J = 1, sNy
       ENDDO                    ! bi = myBxLo, myBxHi
      ENDDO                     ! bj = myByLo, myByHi




#endif /* ALLOW_ADDFLUID */
    RETURN
    END
