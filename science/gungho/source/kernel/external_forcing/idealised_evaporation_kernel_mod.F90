!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Apply an idealised surface evaporative moisture flux.
!> @details For a Cartesian convection-resolving setup, this kernel applies
!>          a prescribed surface heat flux into evaporation as a water
!>          vapour mixing ratio tendency in the lowest model layer.
module idealised_evaporation_kernel_mod

  use argument_mod,    only: arg_type,              &
                             GH_FIELD, GH_REAL,     &
                             GH_READ, GH_READWRITE, &
                             GH_SCALAR,             &
                             CELL_COLUMN
  use constants_mod,   only: i_def, r_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,      only: kernel_type
  use driver_water_constants_mod, only: latent_heat_h2o_condensation

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains metadata for the PSy layer.
  type, public, extends(kernel_type) :: idealised_evaporation_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! dmv_forcing
         arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), & ! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), & ! dz_wtheta
         arg_type(GH_SCALAR, GH_REAL, GH_READ)              & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: idealised_evaporation_code
  end type idealised_evaporation_kernel_type

  public :: idealised_evaporation_code

contains

!> @brief Apply an idealised surface evaporative moisture flux.
!> @param[in]     nlayers           The number of layers
!> @param[in,out] dmv_forcing       Water vapour mixing ratio increment data
!> @param[in]     wetrho_in_wth     Wet density in Wtheta space
!> @param[in]     dz_wtheta         Vertical grid spacing in Wtheta space
!> @param[in]     dt                The model timestep length
!> @param[in]     ndf_wth           Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wth          Number of unique degrees of freedom for Wtheta
!> @param[in]     map_wth           Dofmap for the cell at base of the column for Wtheta
  subroutine idealised_evaporation_code(nlayers,                    &
                                        dmv_forcing,                &
                                        wetrho_in_wth, dz_wtheta,   &
                                        dt,                         &
                                        ndf_wth, undf_wth, map_wth)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dmv_forcing
    real(kind=r_def), dimension(undf_wth), intent(in)    :: wetrho_in_wth, &
                                                            dz_wtheta
    real(kind=r_def),                        intent(in)  :: dt

    integer(kind=i_def), dimension(ndf_wth), intent(in)  :: map_wth

    integer(kind=i_def) :: k

    real(kind=r_def), parameter :: evaporative_heat_flux = 1000.0_r_def ! W m-2

    do k = 0, nlayers
      dmv_forcing(map_wth(1) + k) = 0.0_r_def
    end do

    !> @todo Implement a configurable surface evaporative heat flux, rather
    !>       than the hardwired value used here.
    dmv_forcing(map_wth(1)) = dmv_forcing(map_wth(1)) + (evaporative_heat_flux * dt) / &
        (latent_heat_h2o_condensation * wetrho_in_wth(map_wth(1)) * dz_wtheta(map_wth(1)))

  end subroutine idealised_evaporation_code

end module idealised_evaporation_kernel_mod
