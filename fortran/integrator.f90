MODULE integrator
  USE params, only: ndim
  USE tensor, only:sparse_mul3
  USE aotensor_def, only: aotensor
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1,buf_f0,buf_f1,dW

  INTEGER :: idum1
  
  PUBLIC :: init_integrator, step

CONTAINS
  
  SUBROUTINE init_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim), dW(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    idum1=-1254
    dW=0.D0
  END SUBROUTINE init_integrator
  
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies
  
  SUBROUTINE step(y,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    ! EXTERNAL :: f
    ! INTERFACE
    !    SUBROUTINE f(t,y,res)
    !      REAL(KIND=8), INTENT(IN) :: t
    !      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: y
    !      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: res
    !    END SUBROUTINE f
    ! END INTERFACE
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    
    CALL tendencies(t,y,buf_f0)
    buf_y1 = y+dt*buf_f0
    CALL tendencies(t,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt
    t=t+dt
  END SUBROUTINE step

END MODULE integrator
