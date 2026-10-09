!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Implicit (backward Euler) vertical diffusion of a single column.
!> @details Solves, for boxes i = 1..n,
!>
!>   mass(i) * (x_new(i) - x_old(i)) / dt =
!>       cond(i)   * (x_new(i+1) - x_new(i))
!>     - cond(i-1) * (x_new(i)   - x_new(i-1))
!>
!>          where cond(i) = rho*K/dz is the conductance between boxes i and
!>          i+1, and the flux through the column ends is zero. Every flux
!>          leaves one box and enters its neighbour, so sum(mass*x) is
!>          conserved exactly, and the matrix is diagonally dominant so the
!>          solution is stable for any dt and creates no new extrema.
module vert_smag_tridiag_mod

  use constants_mod, only: i_def, r_def

  implicit none

  private

  public :: implicit_column_diffusion

contains

!> @brief Advance one column of values by one implicit diffusion step.
!> @param[in]     n    Number of boxes in the column
!> @param[in]     dt   Timestep
!> @param[in]     mass Mass per unit area of each box (must be > 0)
!> @param[in]     cond Conductance rho*K/dz between box i and i+1 (>= 0)
!> @param[in,out] x    Values at the start of the step on input, and at the
!>                     end of the step on output
pure subroutine implicit_column_diffusion(n, dt, mass, cond, x)

  implicit none

  integer(kind=i_def), intent(in)    :: n
  real(kind=r_def),    intent(in)    :: dt
  real(kind=r_def),    intent(in)    :: mass(n)
  real(kind=r_def),    intent(in)    :: cond(n-1)
  real(kind=r_def),    intent(inout) :: x(n)

  real(kind=r_def) :: lower(n), diag(n), upper(n)
  real(kind=r_def) :: c_prime(n), d_prime(n)
  real(kind=r_def) :: denom
  integer(kind=i_def) :: i

  if (n < 2) return

  ! Assemble the tridiagonal system, with zero conductance beyond the ends
  do i = 1, n
    lower(i) = 0.0_r_def
    upper(i) = 0.0_r_def
    if (i > 1) lower(i) = -dt * cond(i-1)
    if (i < n) upper(i) = -dt * cond(i)
    diag(i) = mass(i) - lower(i) - upper(i)
  end do

  ! Thomas algorithm: forward elimination ...
  c_prime(1) = upper(1) / diag(1)
  d_prime(1) = mass(1) * x(1) / diag(1)
  do i = 2, n
    denom      = diag(i) - lower(i) * c_prime(i-1)
    c_prime(i) = upper(i) / denom
    d_prime(i) = (mass(i) * x(i) - lower(i) * d_prime(i-1)) / denom
  end do

  ! ... then back substitution
  x(n) = d_prime(n)
  do i = n-1, 1, -1
    x(i) = d_prime(i) - c_prime(i) * x(i+1)
  end do

end subroutine implicit_column_diffusion

end module vert_smag_tridiag_mod
