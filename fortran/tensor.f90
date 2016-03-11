MODULE tensor
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE

  TYPE :: coolist_elem
     INTEGER :: j,k
     REAL(KIND=8) :: v
  END TYPE coolist_elem

  TYPE, PUBLIC :: coolist
     TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: elems
     INTEGER :: nelems = 0
  END TYPE coolist

  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: sparse_mul3,sparse_mul2,copy_tensor,mat_to_coo

CONTAINS
    
  SUBROUTINE copy_tensor(src,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,AllocStat
    
    DO i=1,ndim
       IF (dst(i)%nelems/=0) STOP "*** copy_tensor : Destination coolist not empty ! ***"
       ALLOCATE(dst(i)%elems(src(i)%nelems), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO j=1,src(i)%nelems
          dst(i)%elems(j)%j=src(i)%elems(j)%j
          dst(i)%elems(j)%k=src(i)%elems(j)%k
          dst(i)%elems(j)%v=src(i)%elems(j)%v
       ENDDO
       dst(i)%nelems=src(i)%nelems
    ENDDO
  END SUBROUTINE copy_tensor

  SUBROUTINE mat_to_coo(src,dst)
    REAL(KIND=8), DIMENSION(0:ndim,0:ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,n,AllocStat
    DO i=1,ndim
       n=0
       DO j=1,ndim
          IF (ABS(src(i,j))>real_eps) n=n+1
       ENDDO
       IF (dst(i)%nelems/=0) STOP "*** mat_to_coo : Destination coolist not empty ! ***"
       ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
       DO j=1,ndim
          IF (ABS(src(i,j))>real_eps) THEN
             n=n+1
             dst(i)%elems(n)%j=j
             dst(i)%elems(n)%k=0
             dst(i)%elems(n)%v=src(i,j)
          ENDIF
       ENDDO
       dst(i)%nelems=n
    ENDDO
  END SUBROUTINE mat_to_coo

  SUBROUTINE sparse_mul3(coolist_ijk, arr_j, arr_k, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j, arr_k
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         k=coolist_ijk(i)%elems(n)%k
         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * arr_j(j)*arr_k(k)
      END DO
   END DO
  END SUBROUTINE sparse_mul3

  SUBROUTINE sparse_mul2(coolist_ij, arr_j, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ij(i)%nelems
         j=coolist_ij(i)%elems(n)%j
         res(i) = res(i) + coolist_ij(i)%elems(n)%v * arr_j(j)
      END DO
   END DO
 END SUBROUTINE sparse_mul2
    

END MODULE tensor

