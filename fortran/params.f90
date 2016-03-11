MODULE params

  !-----------------------------------------------------------------------|
  !                                                                       |
  !                  PARAMETERS MODULE                                    |
  !                                                                       |
  !-----------------------------------------------------------------------|---------------|
  !                                                                                       |
  !                                                                                       |
  !      Rem.: Once the init subroutine is called, the parameters are                     |
  !            loaded globally in the main program and its subroutines and functions!     |
  !                                                                                       |
  !---------------------------------------------------------------------------------------|


  IMPLICIT NONE

  PUBLIC

  REAL(KIND=8) :: n,phi0,rra,sig0,k,kp,r,d,f0,gp,H,phi0_npi
  REAL(KIND=8) :: lambda,Co,Go,Ca,To0,Ta0,epsa,Ga,RR

  REAL(KIND=8) :: scale,pi,LR,G,rp,dp,kd,kdp
  REAL(KIND=8) :: Cpo,Lpo,Cpa,Lpa,sBpo,sBpa,LSBpo,LSBpa
  REAL(KIND=8) :: L,sc,sB,sq2
  REAL(KIND=8) :: betp

  REAL(KIND=8) :: t_trans,t_run,dt,tw
  LOGICAL :: writeout

  INTEGER :: nboc,nbatm ! number of blocks
  INTEGER :: natm=0,noc=0 ! number of basis functions
  INTEGER :: ndim ! number of variables
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: oms,ams ! ocean and atmos mode selection

  PRIVATE :: init_nml


  CONTAINS

    
    !-----------------------------------------------------!
    !                                                     !
    ! Read the basic parameters and mode selection from   !
    ! the namelist                                        ! 
    !                                                     !
    !-----------------------------------------------------!

    SUBROUTINE init_nml
      INTEGER :: AllocStat

      NAMELIST /aoscale/  scale,f0,n,rra,phi0_npi
      NAMELIST /oparams/  gp,r,H,d
      NAMELIST /aparams/  k,kp,sig0
      NAMELIST /toparams/ Go,Co,To0
      NAMELIST /taparams/ Ga,Ca,epsa,Ta0
      NAMELIST /otparams/ sc,lambda,RR,sB

      NAMELIST /modeselection/ oms,ams
      NAMELIST /numblocs/ nboc,nbatm

      NAMELIST /int_params/ t_trans,t_run,dt,tw,writeout

      OPEN(8, file="params.nml", status='OLD', recl=80, delim='APOSTROPHE')

      READ(8,nml=aoscale)
      READ(8,nml=oparams)
      READ(8,nml=aparams)
      READ(8,nml=toparams)
      READ(8,nml=taparams)
      READ(8,nml=otparams)
      
      CLOSE(8)

      OPEN(8, file="modeselection.nml", status='OLD', recl=80, delim='APOSTROPHE')
      READ(8,nml=numblocs)

      ALLOCATE(oms(nboc,2),ams(nbatm,2), STAT=AllocStat)
      IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

      READ(8,nml=modeselection)
      CLOSE(8)

      ! Must still sort modeselection here

      OPEN(8, file="int_params.nml", status='OLD', recl=80, delim='APOSTROPHE')
      READ(8,nml=int_params)


    END SUBROUTINE init_nml

    !-----------------------------------------------------!
    !                                                     !
    ! Initialisation routine                              !
    !                                                     !
    !-----------------------------------------------------!

    SUBROUTINE init_params
      INTEGER, DIMENSION(2) :: s
      INTEGER :: i
      CALL init_nml

      !---------------------------------------------------------!
      !                                                         !
      ! Computation of the dimension of the atmospheric         !
      ! and oceanic components                                  !
      !                                                         !
      !---------------------------------------------------------!

      natm=0
      DO i=1,nbatm
         IF (ams(i,1)==1) THEN
            natm=natm+3
         ELSE
            natm=natm+2
         ENDIF
      ENDDO
      s=shape(oms)
      noc=s(1)

      ndim=2*natm+2*noc

      !---------------------------------------------------------!
      !                                                         !
      ! Some general parameters (Domain, beta, gamma, coupling) !
      !                                                         !
      !---------------------------------------------------------!

      pi=dacos(-1.D0)
      L=scale/pi
      phi0=phi0_npi*pi
      LR=sqrt(gp*H)/f0
      G=-L**2/LR**2
      betp=L/rra*cos(phi0)/sin(phi0)
      rp=r/f0
      dp=d/f0
      kd=k*2
      kdp=kp

      !-----------------------------------------------------!
      !                                                     !
      ! DERIVED QUANTITIES                                  !
      !                                                     !
      !-----------------------------------------------------!

      Cpo=Co/(Go*f0) * RR/(f0**2*L**2)
      Lpo=lambda/(Go*f0)
      Cpa=Ca/(Ga*f0) * RR/(f0**2*L**2)/2 ! Cpa acts on psi1-psi3, not on theta
      Lpa=lambda/(Ga*f0)
      sBpo=4*sB*To0**3/(Go*f0) ! long wave radiation lost by ocean to atmosphere space
      sBpa=8*epsa*sB*Ta0**3/(Go*f0) ! long wave radiation from atmosphere absorbed by ocean
      LSBpo=2*epsa*sB*To0**3/(Ga*f0) ! long wave radiation from ocean absorbed by atmosphere
      LSBpa=8*epsa*sB*Ta0**3/(Ga*f0) ! long wave radiation lost by atmosphere to space & ocea

    
    END SUBROUTINE init_params
END MODULE params
