MODULE ic_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr
  USE inprod_analytic, only:awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: IC

  PUBLIC ::load_IC

CONTAINS

  SUBROUTINE load_IC
    INTEGER :: i,AllocStat
    CHARACTER(len=20) :: fm

    NAMELIST /IClist/ IC

    fm(1:6)='(F3.1)'
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(IC(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./IC.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="IC.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=IClist)
       CLOSE(8)
    ELSE
       OPEN(8, file="IC.nml", status='NEW')
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,'(a)') "! Namelist file :                                                              !"
       WRITE(8,'(a)') "! Initial condition.                                                           !"
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,*) ""
       WRITE(8,'(a)') "&ICLIST"
       WRITE(8,*) " ! psi variables"
       DO i=1,natm
          WRITE(8,*) " IC("//TRIM(str(i))//") = 0.D0"//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
       END DO
       WRITE(8,*) " ! theta variables"
       DO i=1,natm
          WRITE(8,*) " IC("//TRIM(str(i+natm))//") = 0.D0"//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
       END DO

       WRITE(8,*) " ! A variables"
       DO i=1,noc
          WRITE(8,*) " IC("//TRIM(str(i+2*natm))//") = 0.D0"//"   ! Nx&
               &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
       END DO
       WRITE(8,*) " ! T variables"
       DO i=1,noc
          WRITE(8,*) " IC("//TRIM(str(i+noc+2*natm))//") = 0.D0"//"   &
               &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
       END DO

       WRITE(8,'(a)') "&END"
       WRITE(8,*) ""
       CLOSE(8)
       WRITE(6,*) "*** IC.nml namelist written. Starting with 0 as initial condition !***"
       IC=0
    ENDIF
    IC(0)=1.0D0
  END SUBROUTINE load_IC
END MODULE ic_def
