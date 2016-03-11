! Generated Fortran90/95 code from aotensor.lua
MODULE aotensor_def

  !---------------------------------------------------------------------------!
  ! aotensor.f90                                                              !
  ! (C) 2013-2014 Lesley De Cruz & Jonathan Demaeyer                          !
  ! See LICENSE.txt for license information.                                  !
  !---------------------------------------------------------------------------!
  !  The equation tensor for the coupled ocean-atmosphere model               !
  !  with temperature which allows for an extensible set of modes             !
  !  in the ocean and in the atmosphere.                                      !
  !---------------------------------------------------------------------------!

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params
  USE inprod_analytic
  USE tensor, only:coolist
  IMPLICIT NONE

  PRIVATE

  INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems

  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: init_aotensor

  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: aotensor


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations                               !
  !                                                     !
  !-----------------------------------------------------!


  FUNCTION psi(i)
    INTEGER :: i,psi
    psi = i
  END FUNCTION psi

  FUNCTION theta(i)
    INTEGER :: i,theta
    theta = i + natm
  END FUNCTION theta

  FUNCTION A(i)
    INTEGER :: i,A
    A = i + 2 * natm
  END FUNCTION A

  FUNCTION T(i)
    INTEGER :: i,T
    T = i + 2 * natm + noc
  END FUNCTION T

  FUNCTION kdelta(i,j)
    INTEGER :: i,j,kdelta
    kdelta=0
    IF (i == j) kdelta = 1
  END FUNCTION kdelta

  SUBROUTINE coeff(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(aotensor)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (.NOT. ALLOCATED(aotensor(i)%elems)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (ABS(v) .ge. real_eps) THEN
       n=(aotensor(i)%nelems)+1
       IF (j .LE. k) THEN
         aotensor(i)%elems(n)%j=j
         aotensor(i)%elems(n)%k=k
       ELSE
         aotensor(i)%elems(n)%j=k
         aotensor(i)%elems(n)%k=j
       END IF
       aotensor(i)%elems(n)%v=v
       aotensor(i)%nelems=n
    END IF
  END SUBROUTINE coeff

  SUBROUTINE simplify_aotensor
    INTEGER :: i,j,k
    INTEGER :: li,lii,liii,n
    IF (.NOT. ALLOCATED(aotensor)) STOP "*** simplify_aotensor routine : tensor not yet allocated ***"
    DO i= 1,ndim
       IF (.NOT. ALLOCATED(aotensor(i)%elems)) STOP "*** simplify_aotensor routine : tensor not yet allocated ***"
       n=aotensor(i)%nelems
       DO li=n,2,-1
          j=aotensor(i)%elems(li)%j
          k=aotensor(i)%elems(li)%k
          DO lii=li-1,1,-1
             IF ((j==aotensor(i)%elems(lii)%j).AND.(k==aotensor(i)%elems(lii)%k)) THEN
                ! Found another entry with the same i,j,k: merge both into
                ! the one listed first (of those two). 
                aotensor(i)%elems(lii)%v=aotensor(i)%elems(lii)%v+aotensor(i)%elems(li)%v
                ! Shift the rest of the items one place down.
                DO liii=li+1,n
                   aotensor(i)%elems(liii-1)%j=aotensor(i)%elems(liii)%j
                   aotensor(i)%elems(liii-1)%k=aotensor(i)%elems(liii)%k
                   aotensor(i)%elems(liii-1)%v=aotensor(i)%elems(liii)%v
                END DO
                aotensor(i)%nelems=aotensor(i)%nelems-1
                ! Here we should stop because the li no longer points to the
                ! original i,j,k element
                EXIT
             ENDIF
          ENDDO
       ENDDO
       n=aotensor(i)%nelems
       DO li=1,n
         ! Clear new "almost" zero entries and shift rest of the items one place down.
         ! Make sure not to skip any entries while shifting!
          DO WHILE (ABS(aotensor(i)%elems(li)%v) < real_eps)
             DO liii=li+1,n
                aotensor(i)%elems(liii-1)%j=aotensor(i)%elems(liii)%j
                aotensor(i)%elems(liii-1)%k=aotensor(i)%elems(liii)%k
                aotensor(i)%elems(liii-1)%v=aotensor(i)%elems(liii)%v
             ENDDO
             aotensor(i)%nelems=aotensor(i)%nelems-1
          ENDDO
       ENDDO

    ENDDO
  END SUBROUTINE simplify_aotensor
                
  SUBROUTINE add_count(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (ABS(v) .ge. real_eps) count_elems(i)=count_elems(i)+1
  END SUBROUTINE add_count

  SUBROUTINE compute_aotensor(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k
    CALL func(theta(1),0,0,(Cpa / (1 - atmos%a(1,1) * sig0)))
    DO i = 1, natm
       DO j = 1, natm
          CALL func(psi(i),psi(j),0,-(((atmos%c(i,j) * betp) / atmos%a(i,i))) - (kd * kdelta(i,j)) / 2)
          CALL func(theta(i),psi(j),0,(atmos%a(i,j) * kd * sig0) / (-2 + 2 * atmos%a(i,i) * sig0))
          CALL func(psi(i),theta(j),0,(kd * kdelta(i,j)) / 2)
          CALL func(theta(i),theta(j),0,(-((sig0 * (2. * atmos%c(i,j) * betp +&
               & atmos%a(i,j) * (kd + 4. * kdp)))) + 2. * (LSBpa + sc * Lpa) &
               &* kdelta(i,j)) / (-2. + 2. * atmos%a(i,i) * sig0))
          DO k = 1, natm
             CALL func(psi(i),psi(j),psi(k),-((atmos%b(i,j,k) / atmos%a(i,i))))
             CALL func(psi(i),theta(j),theta(k),-((atmos%b(i,j,k) / atmos%a(i,i))))
             CALL func(theta(i),psi(j),theta(k),(atmos%g(i,j,k) -&
                  & atmos%b(i,j,k) * sig0) / (-1 + atmos%a(i,i) *&
                  & sig0))
             CALL func(theta(i),theta(j),psi(k),(atmos%b(i,j,k) * sig0) / (1 - atmos%a(i,i) * sig0))
          END DO
       END DO
       DO j = 1, noc
          CALL func(psi(i),A(j),0,kd * atmos%d(i,j) / (2 * atmos%a(i,i)))
          CALL func(theta(i),A(j),0,kd * (atmos%d(i,j) * sig0) / (2 - 2 * atmos%a(i,i) * sig0))
          CALL func(theta(i),T(j),0,atmos%s(i,j) * (2 * LSBpo + Lpa) / (2 - 2 * atmos%a(i,i) * sig0))
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, natm
          CALL func(A(i),psi(j),0,ocean%K(i,j) * dp / (ocean%M(i,i) + G))
          CALL func(A(i),theta(j),0,-(ocean%K(i,j)) * dp / (ocean%M(i,i) + G))
       END DO
       DO j = 1, noc
          CALL func(A(i),A(j),0,-((ocean%N(i,j) * betp + ocean%M(i&
               &,i) * (rp + dp) * kdelta(i,j))) / (ocean%M(i,i) + G))
          DO k = 1, noc
             CALL func(A(i),A(j),A(k),-(ocean%C(i,j,k)) / (ocean%M(i,i) + G))
          END DO
       END DO
    END DO
    DO i = 1, noc
       CALL func(T(i),0,0,Cpo * ocean%W(i,1))
       DO j = 1, natm
          CALL func(T(i),theta(j),0,ocean%W(i,j) * (2 * sc * Lpo + sBpa))
       END DO
       DO j = 1, noc
          CALL func(T(i),T(j),0,-((Lpo + sBpo)) * kdelta(i,j))
          DO k = 1, noc
             CALL func(T(i),A(j),T(k),-(ocean%O(i,j,k)))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_aotensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  SUBROUTINE init_aotensor
    INTEGER :: i
    INTEGER :: AllocStat 

    CALL init_params  ! Iniatialise the parameter

    CALL init_inprod  ! Initialise the inner product tensors

    ALLOCATE(aotensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0

    CALL compute_aotensor(add_count)

    DO i=1,ndim
       ALLOCATE(aotensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_aotensor(coeff)

    CALL simplify_aotensor

    CALL deallocate_inprod ! Clean the inner product tensors

  END SUBROUTINE init_aotensor
END MODULE aotensor_def
      


