MODULE util
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: str,rstr

CONTAINS
    
  ! SUBROUTINE scalar_allocate(x)
  !   INTEGER :: AllocStat
  !   IF (.NOT. ALLOCATED(x)) THEN 
  !      ALLOCATE(x, STAT=AllocStat)
  
  CHARACTER(len=20) FUNCTION str(k) !   Convert an integer to string.
    INTEGER, INTENT(IN) :: k
    WRITE (str, *) k
    str = ADJUSTL(str)
  END FUNCTION str

  CHARACTER(len=40) FUNCTION rstr(x,fm) !   Convert an integer to string.
    REAL(KIND=8), INTENT(IN) :: x
    CHARACTER(len=20), INTENT(IN) :: fm
    WRITE (rstr, TRIM(ADJUSTL(fm))) x
    rstr = ADJUSTL(rstr)
  END FUNCTION rstr   
    

END MODULE util
