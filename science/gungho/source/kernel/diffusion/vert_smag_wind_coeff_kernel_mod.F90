!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Map the vertical Smagorinsky momentum coefficients onto cell faces.
!> @details Horizontal winds live on the W2 horizontal dofs at layer heights,
!>          so their vertical fluxes are on theta levels at the faces. Each
!>          face takes the average of the two cells either side (weighted by
!>          w2_rmultiplicity) of:
!>            - the conductance rho_wet*K_m/dz at theta level k, stored in
!>              W2 index k for k = 1..nlayers-1 (index 0, the ground, is left
!>              at zero, giving zero flux there);
!>            - the wet mass per unit area of layer k, stored in W2 index k.
!>          Both output fields must be zeroed before this kernel is called.
module vert_smag_wind_coeff_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_INC,   &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: vert_smag_wind_coeff_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),      & ! rhokm_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),      & ! mass_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),  & ! visc_m
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),  & ! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),      & ! wetrho_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),  & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),      & ! height_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2)       & ! w2_rmultiplicity
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vert_smag_wind_coeff_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vert_smag_wind_coeff_code

contains

!> @brief Add this cell's share of the face coefficients.
!> @param[in]     nlayers          Number of layers in the mesh
!> @param[in,out] rhokm_w2         Conductance rho_wet*K_m/dz on faces
!> @param[in,out] mass_w2          Wet mass per unit area of each layer on faces
!> @param[in]     visc_m           Momentum eddy diffusivity K_m on theta levels
!> @param[in]     wetrho_in_wth    Wet density on theta levels
!> @param[in]     wetrho_in_w3     Wet density on layers
!> @param[in]     height_wth       Height of theta levels above the surface
!> @param[in]     height_w3        Height of layer centres above the surface
!> @param[in]     w2_rmultiplicity Reciprocal of the number of cells sharing
!>                                 each W2 dof
!> @param[in]     ndf_w2           Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2          Number of unique degrees of freedom for W2
!> @param[in]     map_w2           Dofmap for the cell at the base of the column
!> @param[in]     ndf_wt           Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wt          Number of unique degrees of freedom for Wtheta
!> @param[in]     map_wt           Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3           Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3          Number of unique degrees of freedom for W3
!> @param[in]     map_w3           Dofmap for the cell at the base of the column
subroutine vert_smag_wind_coeff_code( nlayers,                 &
                                      rhokm_w2, mass_w2,       &
                                      visc_m,                  &
                                      wetrho_in_wth,           &
                                      wetrho_in_w3,            &
                                      height_wth, height_w3,   &
                                      w2_rmultiplicity,        &
                                      ndf_w2, undf_w2, map_w2, &
                                      ndf_wt, undf_wt, map_wt, &
                                      ndf_w3, undf_w3, map_w3  &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(undf_w2), intent(inout) :: rhokm_w2, mass_w2
  real(kind=r_def), dimension(undf_wt), intent(in)    :: visc_m
  real(kind=r_def), dimension(undf_wt), intent(in)    :: wetrho_in_wth
  real(kind=r_def), dimension(undf_w3), intent(in)    :: wetrho_in_w3
  real(kind=r_def), dimension(undf_wt), intent(in)    :: height_wth
  real(kind=r_def), dimension(undf_w3), intent(in)    :: height_w3
  real(kind=r_def), dimension(undf_w2), intent(in)    :: w2_rmultiplicity

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def)    :: layer_mass(0:nlayers-1), conductance(nlayers-1)

  ! Wet mass per unit area of each layer
  do k = 0, nlayers - 1
    layer_mass(k) = wetrho_in_w3(map_w3(1) + k)                          &
                  * (height_wth(map_wt(1) + k + 1) - height_wth(map_wt(1) + k))
  end do

  ! Conductance across theta level k, between layers k-1 and k
  do k = 1, nlayers - 1
    conductance(k) = wetrho_in_wth(map_wt(1) + k) * visc_m(map_wt(1) + k) &
                   / (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + k - 1))
  end do

  ! Horizontal faces only
  do df = 1, 4
    do k = 0, nlayers - 1
      mass_w2(map_w2(df) + k) = mass_w2(map_w2(df) + k)                  &
                              + w2_rmultiplicity(map_w2(df) + k) * layer_mass(k)
    end do
    do k = 1, nlayers - 1
      rhokm_w2(map_w2(df) + k) = rhokm_w2(map_w2(df) + k)                &
                               + w2_rmultiplicity(map_w2(df) + k) * conductance(k)
    end do
  end do

end subroutine vert_smag_wind_coeff_code

end module vert_smag_wind_coeff_kernel_mod
