
! test_tl_ad.f90
!
!> Tests for the Tangent Linear (TL) and Adjoint (AD) model versions
!> of MAOOAM.
!     
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!


PROGRAM test_tl_ad
  USE params, only:ndim,dt,t_trans
  USE aotensor_def, only: init_aotensor
  USE integrator, only: init_integrator,step
  USE tl_ad_tensor, only: init_tltensor, init_adtensor
  USE tl_ad_integrator, only: init_tl_ad_integrator,tl_step,ad_step
  USE stoch_mod, only: gasdev
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y0_IC,y0,y0prime,dy0,dy0_bis
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y1,y1prime
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dy,dy_bis
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dy1,dy1_tl,dy1_bis_tl,dy1_ad,dy1_bis_ad
  REAL(KIND=8) :: t=0.D0
  REAL(KIND=8) :: norm1,norm2
  INTEGER :: i,n
  
  ! Compute the tensors

  CALL init_aotensor
  CALL init_tltensor
  CALL init_adtensor

  CALL init_integrator  ! Initialize the model integrator
  CALL init_tl_ad_integrator  ! Initialize the TL & AD integrator

  ALLOCATE(y0_IC(0:ndim),dy(0:ndim),y0(0:ndim),y0prime(0:ndim)&
       &,y1(0:ndim))!, STAT=AllocStat)
  ALLOCATE(y1prime(0:ndim),dy0(0:ndim),dy0_bis(0:ndim))!, STAT&
       ! &=AllocStat)
  ALLOCATE(dy1(0:ndim),dy1_tl(0:ndim),dy_bis(0:ndim)&
       &,dy1_bis_tl(0:ndim),dy1_ad(0:ndim),dy1_bis_ad(0:ndim))!, STAT&
       ! &=AllocStat)
 

  ! Test Taylor property for the Tangent Linear.
  ! lim(\lambda->0) M(x+\lambda dx) - M(x) / M'(\lambda dx) = 1

  y0_IC(0)=1.D0

  ! Set all values to random.

  DO i=1,ndim
     y0_IC(i)=0.01*gasdev()
  ENDDO

  ! Evolve during transient period

  ! PRINT*, 'Random values:',y0_IC(0:ndim)
  DO WHILE (t<t_trans)
     CALL step(y0_IC,t,dt,y0)
     y0_IC=y0
  END DO
  PRINT*, 'Initial values:',y0_IC(0:ndim)

  ! Test 1: Taylor test
  ! Integrate the original model by one step, integrate with perturbed
  ! IC, and test if the difference approximates the TL

  DO n=0,6
     ! Small perturbation.
     dy=2.D0**(-n)/sqrt(real(ndim))
     dy(0)=0.D0
     PRINT*, "Perturbation size:",dot_product(dy,dy)
     
     y0 = y0_IC*1
     y0prime = y0 + dy
     CALL step(y0,t,dt,y1)
     CALL step(y0prime,t,dt,y1prime)

     dy1 = y1prime - y1

     dy0 = dy*1
     CALL tl_step(dy0,y0_IC,t,dt,dy1_tl)

     ! Don't forget to set 0'th component to 0...
     dy1(0)=0.D0
     dy1_tl(0)=0.D0

     PRINT*, "Resulting difference in trajectory: (epsilon ~ 2^-",n,")"
     PRINT*, "diff:    ",dot_product(dy1,dy1)
     PRINT*, "tl:      ",dot_product(dy1_tl,dy1_tl)
     PRINT*, "ratio:   ",dot_product(dy1,dy1)/dot_product(dy1_tl,dy1_tl)
  END DO

  ! Test 2: Adjoint Identity: <M(TL).x,y> = <x,M(AD).y>
  
  DO i=1,100
     ! Any perturbation.
     DO n=1,ndim
        dy(n)=gasdev()
     END DO
     dy(0)=0.D0

     DO n=1,ndim
        dy_bis(n)=gasdev()
     END DO
     dy_bis(0)=0.D0

     ! Calculate M(TL).x in dy1_tl
     dy0 = dy*1
     CALL tl_step(dy0,y0_IC,t,dt,dy1_tl)
     
     ! Calculate M(AD).x in dy1_ad
     CALL ad_step(dy0,y0_IC,t,dt,dy1_ad)

     ! Calculate M(TL).y in dy1_bis_tl
     CALL tl_step(dy0_bis,y0_IC,t,dt,dy1_bis_tl)
     
     ! Calculate M(AD).y in dy1_bis_ad
     CALL ad_step(dy0_bis,y0_IC,t,dt,dy1_bis_ad)

     ! Calculate norm <M(TL).x,y>
     norm1 = dot_product(dy1_tl,dy0_bis)
     ! Calculate norm <x,M(AD).y>
     norm2 = dot_product(dy0,dy1_bis_ad)

     PRINT*, "<M(TL).x,y> = ", norm1
     PRINT*, "<x,M(AD).y> = ", norm2
     PRINT*, "Ratio       = ", norm1/norm2

     ! Calculate norm <M(TL).y,x>
     norm1 = dot_product(dy1_bis_tl,dy0)
     ! Calculate norm <y,M(AD).x>
     norm2 = dot_product(dy0_bis,dy1_ad)

     PRINT*, "<M(TL).y,x> = ", norm1
     PRINT*, "<y,M(AD).x> = ", norm2
     PRINT*, "Ratio       = ", norm1/norm2
  END DO

END PROGRAM test_tl_ad

