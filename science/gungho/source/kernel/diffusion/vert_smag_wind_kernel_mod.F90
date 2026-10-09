!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Implicit vertical Smagorinsky diffusion of the horizontal wind.
!> @details For each horizontal face of the cell, solves
!>          d(u)/dt = (1/rho_wet) d/dz (rho_wet K_m du/dz) with backward Euler
!>          over the layers, with zero flux at the ground and the model top.
!>          The face coefficients come from vert_smag_wind_coeff_kernel_mod.
!>          Every value written on a face depends only on that face, so the
!>          two cells sharing it write identical increments. Vertical wind
!>          dofs are not mixed and get a zero increment.
module vert_smag_wind_kernel_mod

  use argument_mod,          only : arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_SCALAR,         &
                                    GH_READ, GH_WRITE, &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2
  use kernel_mod,            only : kernel_type
  use vert_smag_tridiag_mod, only : implicit_column_diffusion

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: vert_smag_wind_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/              &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2), & ! du
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2), & ! u_phys
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2), & ! rhokm_w2
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2), & ! mass_w2
         arg_type(GH_SCALAR, GH_REAL, GH_READ)       & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vert_smag_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vert_smag_wind_code

contains

!> @brief Compute the wind increment from one implicit vertical diffusion step.
!> @param[in]     nlayers  Number of layers in the mesh
!> @param[in,out] du       Wind increment (m/s, not scaled by face area)
!> @param[in]     u_phys   Wind normal to each face (m/s)
!> @param[in]     rhokm_w2 Conductance rho_wet*K_m/dz below layer k on faces
!> @param[in]     mass_w2  Wet mass per unit area of layer k on faces
!> @param[in]     dt       The model timestep length
!> @param[in]     ndf_w2   Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2  Number of unique degrees of freedom for W2
!> @param[in]     map_w2   Dofmap for the cell at the base of the column
subroutine vert_smag_wind_code( nlayers,                &
                                du, u_phys,             &
                                rhokm_w2, mass_w2,      &
                                dt,                     &
                                ndf_w2, undf_w2, map_w2 &
                              )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_w2), intent(inout) :: du
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u_phys
  real(kind=r_def), dimension(undf_w2), intent(in)    :: rhokm_w2, mass_w2
  real(kind=r_def),                     intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def)    :: mass(nlayers), cond(nlayers-1), x(nlayers)

  ! Horizontal faces: box i of the solver is layer i-1, and cond(i) is the
  ! conductance across theta level i, between layers i-1 and i
  do df = 1, 4
    do k = 1, nlayers
      mass(k) = mass_w2(map_w2(df) + k - 1)
      x(k)    = u_phys(map_w2(df) + k - 1)
    end do
    do k = 1, nlayers - 1
      cond(k) = rhokm_w2(map_w2(df) + k)
    end do

    call implicit_column_diffusion(nlayers, dt, mass, cond, x)

    do k = 1, nlayers
      du(map_w2(df) + k - 1) = x(k) - u_phys(map_w2(df) + k - 1)
    end do
  end do

  ! Vertical faces are not mixed
  do df = 5, 6
    do k = 0, nlayers - 1
      du(map_w2(df) + k) = 0.0_r_def
    end do
  end do

end subroutine vert_smag_wind_code

end module vert_smag_wind_kernel_mod
