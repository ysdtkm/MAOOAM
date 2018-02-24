PROGRAM run_true_and_tlm
  USE params, only:ndim,dt,t_trans,t_run,tw,writeout
  USE IC_def, only: load_IC, IC
  USE aotensor_def, only: init_aotensor
  USE integrator, only: init_integrator,step
  USE tl_ad_tensor, only: init_tltensor, init_adtensor
  USE tl_ad_integrator, only: init_tl_ad_integrator,tl_step,ad_step
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y0_IC,y0
  REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: tlm
  REAL(KIND=8) :: t=0.D0
  INTEGER :: irec, i, j

  CALL init_aotensor
  CALL load_IC
  CALL init_tltensor
  CALL init_adtensor
  CALL init_integrator
  CALL init_tl_ad_integrator

  ALLOCATE (y0_IC(0:ndim), tlm(0:ndim, 0:ndim))

  y0_IC = IC

  ! Evolve during transient period
  DO WHILE (t < t_trans)
     CALL step(y0_IC, t, dt, y0)
     y0_IC = y0
  END DO


  IF (writeout) OPEN(10,file='evol_field_tlm.dat')
  DO WHILE (t < t_run)
    print *, t
    CALL tl_matrix_analytic(y0_IC, t, dt, int(tw / dt), tlm)
    IF (writeout) WRITE(10,*) t, y0_IC(1:ndim), tlm(1:ndim, 1:ndim)
  END DO
  IF (writeout) CLOSE(10)

!  OPEN(10, file='evol_field_tlm.dat', form="unformatted", access="direct", recl=8)
!  irec = 1
!  DO WHILE (t < t_run)
!    print *, t
!    CALL tl_matrix_analytic(y0_IC, t, dt, int(tw / dt), tlm)
!    DO i = 1, ndim
!      WRITE(10, rec=irec) y0_IC(i); irec = irec + 1
!    END DO
!    DO i = 1, ndim
!      DO j = 1, ndim
!        WRITE(10, rec=irec) tlm(j, i); irec = irec + 1
!      END DO
!    END DO
!    ! WRITE(10,*) t, y0_IC(1:ndim), tlm(1:ndim, 1:ndim)
!  END DO
!  CLOSE(10)

  DEALLOCATE (y0_IC, tlm)

END PROGRAM run_true_and_tlm

SUBROUTINE tl_matrix_analytic(xctl, t, dt, nt, tlm)
  USE params, only: ndim
  USE integrator, only: step
  USE tl_ad_integrator, only: tl_step
  IMPLICIT NONE
  REAL(8), INTENT(INOUT) :: xctl(0:ndim), t
  REAL(8), INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: nt
  REAL(8), INTENT(OUT) :: tlm(0:ndim, 0:ndim)

  INTEGER :: i, j, k
  REAL(8) :: t_dummy = 0.0, tmp(0:ndim)

  tlm(:, :) = 0.0D0
  do i = 0, ndim
    tlm(i, i) = 1.0D0
  end do

  do k = 1, nt
    do j = 0, ndim
      CALL tl_step(tlm(:, j), xctl, t_dummy, dt, tmp(:))
      tlm(:, j) = tmp(:)
    end do
    CALL step(xctl, t, dt, tmp)
    xctl(:) = tmp(:)
  end do
END SUBROUTINE tl_matrix_analytic


