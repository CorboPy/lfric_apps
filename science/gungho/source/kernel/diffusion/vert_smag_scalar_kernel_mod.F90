!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Implicit vertical Smagorinsky diffusion of a scalar held in Wtheta.
!> @details Solves d(q)/dt = (1/rho) d/dz (rho K_h dq/dz) with backward Euler
!>          on theta levels 1..nlayers, with zero flux at the ground and the
!>          model top. The ground value (level 0) is not changed.
!>
!>          Level k is weighted by the dry mass 0.5*(rho(k-1)*dz(k-1) +
!>          rho(k)*dz(k)), the same weighting moisture_conservation_alg uses,
!>          so the column total it reports is conserved exactly. This assumes
!>          the column cross-section does not vary with height (planar or
!>          shallow-atmosphere meshes).
module vert_smag_scalar_kernel_mod

  use argument_mod,          only : arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_SCALAR,         &
                                    GH_READ, GH_WRITE, &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W3
  use kernel_mod,            only : kernel_type
  use vert_smag_tridiag_mod, only : implicit_column_diffusion

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: vert_smag_scalar_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                   &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),  & ! increment
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),  & ! field
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),  & ! visc_h
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),      & ! rho
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),  & ! height_wth
         arg_type(GH_SCALAR, GH_REAL, GH_READ)            & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vert_smag_scalar_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vert_smag_scalar_code

contains

!> @brief Compute the increment from one implicit vertical diffusion step.
!> @param[in]     nlayers    Number of layers in the mesh
!> @param[in,out] increment  Increment to the scalar over the step
!> @param[in]     field      Scalar to be diffused
!> @param[in]     visc_h     Scalar eddy diffusivity K_h on theta levels
!> @param[in]     rho        Dry density
!> @param[in]     height_wth Height of theta levels above the surface
!> @param[in]     dt         The model timestep length
!> @param[in]     ndf_wt     Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wt    Number of unique degrees of freedom for Wtheta
!> @param[in]     map_wt     Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3    Number of unique degrees of freedom for W3
!> @param[in]     map_w3     Dofmap for the cell at the base of the column
subroutine vert_smag_scalar_code( nlayers,                 &
                                  increment,               &
                                  field,                   &
                                  visc_h,                  &
                                  rho,                     &
                                  height_wth,              &
                                  dt,                      &
                                  ndf_wt, undf_wt, map_wt, &
                                  ndf_w3, undf_w3, map_w3  &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(undf_wt), intent(inout) :: increment
  real(kind=r_def), dimension(undf_wt), intent(in)    :: field
  real(kind=r_def), dimension(undf_wt), intent(in)    :: visc_h
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_wt), intent(in)    :: height_wth
  real(kind=r_def),                     intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: dz(0:nlayers-1), rho_dz(0:nlayers-1)
  real(kind=r_def)    :: mass(nlayers), cond(nlayers-1), x(nlayers)

  ! Layer thicknesses and dry mass per unit area of each layer
  do k = 0, nlayers - 1
    dz(k)     = height_wth(map_wt(1) + k + 1) - height_wth(map_wt(1) + k)
    rho_dz(k) = rho(map_w3(1) + k) * dz(k)
  end do

  ! Unknowns are theta levels 1..nlayers. Each level owns half of the layer
  ! below and half of the layer above; the top level only has a layer below.
  do k = 1, nlayers - 1
    mass(k) = 0.5_r_def * (rho_dz(k-1) + rho_dz(k))
  end do
  mass(nlayers) = 0.5_r_def * rho_dz(nlayers-1)

  ! Conductance rho*K_h/dz across layer k, between theta levels k and k+1.
  ! There is no conductance across layer 0, giving zero flux at the ground.
  do k = 1, nlayers - 1
    cond(k) = rho(map_w3(1) + k) * 0.5_r_def                          &
            * (visc_h(map_wt(1) + k) + visc_h(map_wt(1) + k + 1)) / dz(k)
  end do

  do k = 1, nlayers
    x(k) = field(map_wt(1) + k)
  end do

  call implicit_column_diffusion(nlayers, dt, mass, cond, x)

  increment(map_wt(1)) = 0.0_r_def
  do k = 1, nlayers
    increment(map_wt(1) + k) = x(k) - field(map_wt(1) + k)
  end do

end subroutine vert_smag_scalar_code

end module vert_smag_scalar_kernel_mod
