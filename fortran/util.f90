

! util.f90
!
!>  Utility module
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE util
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: str,rstr

CONTAINS
    
  ! SUBROUTINE scalar_allocate(x)
  !   INTEGER :: AllocStat
  !   IF (.NOT. ALLOCATED(x)) THEN 
  !      ALLOCATE(x, STAT=AllocStat)
  
  !>  Convert an integer to string.
  CHARACTER(len=20) FUNCTION str(k) 
    INTEGER, INTENT(IN) :: k
    WRITE (str, *) k
    str = ADJUSTL(str)
  END FUNCTION str

  !>   Convert a real to string with a given format
  CHARACTER(len=40) FUNCTION rstr(x,fm) 
    REAL(KIND=8), INTENT(IN) :: x
    CHARACTER(len=20), INTENT(IN) :: fm
    WRITE (rstr, TRIM(ADJUSTL(fm))) x
    rstr = ADJUSTL(rstr)
  END FUNCTION rstr   
    

END MODULE util
