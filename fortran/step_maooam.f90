subroutine step_maooam(x0, dt)
  use params, only: ndim
  use aotensor_def, only: init_aotensor
  use integrator, only: init_integrator,step

  implicit none

  real(8), intent(inout) :: x0(0:ndim)
  real(8), intent(in) :: dt
  real(8) :: t_dummy
  real(8), save, allocatable :: x1(:)
  logical, save :: first_time = .true.
  integer :: stat

  if (first_time) then
    call init_aotensor
    call init_integrator
    first_time = .false.
    t_dummy = 0.0d0
    allocate (x1(0:ndim), stat=stat)
    if (stat /= 0) stop "*** Allocation error of x1 ! ***"
  end if

  x1(:) = 0.0
  call step(x0, t_dummy, dt, x1)
  x0 = x1
end subroutine step_maooam

