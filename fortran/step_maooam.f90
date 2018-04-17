subroutine step_maooam(x0, dt)
  use params, only: ndim
  use aotensor_def, only: init_aotensor
  use integrator, only: init_integrator,step

  implicit none

  real(8), intent(inout) :: x0(ndim)
  real(8), intent(in) :: dt
  real(8) :: t_dummy
  logical, save :: first_time = .true.

  call init_aotensor
  call init_integrator
  first_time = .false.
  t_dummy = 0.0d0

  call step(x0, t_dummy, dt, x0)
end subroutine step_maooam

