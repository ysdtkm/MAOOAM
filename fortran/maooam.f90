!  maooam.f90
!  (C) 2015 Lesley De Cruz & Jonathan Demaeyer
!  See LICENSE.txt for license information.

!---------------------------------------------------------------------------!
! Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere !
! model MAOOAM.                                                             !
!                                                                           !
!---------------------------------------------------------------------------!


PROGRAM maooam 
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X,Xnew
  REAL(KIND=8) :: t=0.D0

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  IF (writeout) OPEN(10,file='evol_field.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  X=IC

  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0

  DO WHILE (t<t_run)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        CALL acc(X)
     END IF
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam 
