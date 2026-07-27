!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Apply idealised tropospheric cooling and stratospheric nudging.
!> @details For a Cartesian convection-resolving setup, this kernel applies:
!>          1) A column-integrated cooling flux distributed as a uniform
!>             temperature tendency from the surface to a diagnosed tropopause.
!>          2) Newtonian nudging above that level towards a fixed tropopause
!>             temperature on a fixed timescale.
module idealised_radiation_kernel_mod

  use argument_mod,    only: arg_type,              &
                             GH_FIELD, GH_REAL,     &
                             GH_READ, GH_READWRITE, &
                             GH_SCALAR,             &
                             CELL_COLUMN
  use constants_mod,   only: i_def, r_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,      only: kernel_type
  use planet_config_mod, only: cp, gravity, p_zero, kappa

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains metadata for the PSy layer.
  type, public, extends(kernel_type) :: idealised_radiation_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                   &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: idealised_radiation_code
  end type idealised_radiation_kernel_type

  public :: idealised_radiation_code

contains

!> @brief Apply idealised tropospheric cooling and stratospheric nudging.
!> @param[in]     nlayers The number of layers
!> @param[in,out] dtheta Potential temperature increment data
!> @param[in]     theta Potential temperature data
!> @param[in]     exner_in_wth Exner pressure in Wtheta space
!> @param[in]     temperature_mean Horizontally averaged absolute temperature
!> @param[in]     dt The model timestep length
!> @param[in]     ndf_wth Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wth Number of unique degrees of freedom for Wtheta
!> @param[in]     map_wth Dofmap for the cell at base of the column for Wtheta
  subroutine idealised_radiation_code(nlayers,                    &
                                      dtheta, theta,              &
                                      exner_in_wth,               &
                                      temperature_mean, dt,       &
                                      ndf_wth, undf_wth, map_wth  &
                                     )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta
    real(kind=r_def), dimension(undf_wth), intent(in)    :: theta
    real(kind=r_def), dimension(undf_wth), intent(in)    :: exner_in_wth
    real(kind=r_def), dimension(undf_wth), intent(in)    :: temperature_mean
    real(kind=r_def),                      intent(in)     :: dt

    integer(kind=i_def), dimension(ndf_wth), intent(in)  :: map_wth

    integer(kind=i_def) :: k, tropopause_level
    real(kind=r_def)    :: exner, temperature, dtemp_dt
    real(kind=r_def)    :: p_surface, p_tropopause, delta_pressure
    real(kind=r_def)    :: tropospheric_dtemp_dt

    real(kind=r_def), parameter :: column_cooling_flux = 200.0_r_def  ! W m-2
    real(kind=r_def), parameter :: tropopause_temperature = 200.0_r_def ! K
    real(kind=r_def), parameter :: nudging_timescale = 21600.0_r_def ! 6 hours

    tropopause_level = nlayers
    do k = 0, nlayers
      if (temperature_mean(map_wth(1) + k) <= tropopause_temperature) then
        tropopause_level = k
        exit
      end if
    end do

    p_surface = p_zero * exner_in_wth(map_wth(1)) ** (1.0_r_def / kappa)
    p_tropopause = p_zero * exner_in_wth(map_wth(1) + tropopause_level) **        &
                   (1.0_r_def / kappa)
    delta_pressure = p_surface - p_tropopause

    if (delta_pressure > epsilon(1.0_r_def)) then
      tropospheric_dtemp_dt = -column_cooling_flux * gravity / (cp * delta_pressure)
    else
      tropospheric_dtemp_dt = 0.0_r_def
    end if

    do k = 0, nlayers
      exner = exner_in_wth(map_wth(1) + k)

      if (exner <= epsilon(1.0_r_def)) then
        dtheta(map_wth(1) + k) = 0.0_r_def
        cycle
      end if

      if (k <= tropopause_level) then
        dtemp_dt = tropospheric_dtemp_dt
      else
        temperature = theta(map_wth(1) + k) * exner
        dtemp_dt = -(temperature - tropopause_temperature) / nudging_timescale
      end if

      dtheta(map_wth(1) + k) = dtemp_dt * dt / exner
    end do

  end subroutine idealised_radiation_code

end module idealised_radiation_kernel_mod
