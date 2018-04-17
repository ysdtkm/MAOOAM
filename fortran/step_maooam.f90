subroutine step_maooam(x0, dt)
  use params, only: ndim
  use aotensor_def, only: init_aotensor
  use integrator, only: init_integrator,step

  implicit none

  real(8), intent(inout) :: x0(:)
  real(8), intent(in) :: dt
  real(8), allocatable :: x1(:)
  real(8) :: t_dummy
  logical, save :: first_time = .true.

  allocate(x1(ndim))

  if (first_time) then
    call init_aotensor
    call init_integrator
    first_time = .false.
    t_dummy = 0.0d0
  end if

  call step(x0, t_dummy, dt, x1)
  x0 = x1

  deallocate(x1)
end subroutine step_maooam

