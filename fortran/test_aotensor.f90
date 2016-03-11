  !---------------------------------------------------------------------------!
  ! inprod_analytic_test.f90                                                  !
  ! (C) 2013-2014 Lesley De Cruz & Jonathan Demaeyer                          !
  ! See LICENSE.txt for license information.                                  !
  !---------------------------------------------------------------------------!
  !  Small progrem to print the aotensor.                                     !
  !---------------------------------------------------------------------------!

PROGRAM test_aotensor

  USE params, only: ndim
  USE aotensor_def, only: aotensor, init_aotensor
  USE util, only: str
  
  IMPLICIT NONE
  INTEGER :: i,j,k,n 
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16
  
  ! Program

  CALL init_aotensor    ! Compute the tensor

  DO i=1,ndim
    DO n=1,aotensor(i)%nelems
      j=aotensor(i)%elems(n)%j
      k=aotensor(i)%elems(n)%k
      IF( ABS(aotensor(i)%elems(n)%v) .GE. real_eps) THEN
        write(*,"(A,ES12.5)") "aotensor["//TRIM(str(i))//"]["//TRIM(str(j)) &
        &//"]["//TRIM(str(k))//"] = ",aotensor(i)%elems(n)%v
      END IF
    END DO
  END DO

END PROGRAM test_aotensor
