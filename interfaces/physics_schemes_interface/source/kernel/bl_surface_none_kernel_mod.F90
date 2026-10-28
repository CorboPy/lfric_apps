!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Lower boundary values for the UM boundary layer when no surface
!>        scheme is used.
!> @details Sets the surface quantities that the explicit JULES scheme would
!>          otherwise pass to the UM boundary layer, so that the boundary layer
!>          sees a free-slip, zero-flux lower boundary. Only the surface (level
!>          0) values of the boundary layer fields are written.

module bl_surface_none_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_INTEGER,                &
                                GH_READ, GH_WRITE,         &
                                GH_READWRITE, CELL_COLUMN, &
                                ANY_DISCONTINUOUS_SPACE_1, &
                                ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,     only : i_def, r_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: bl_surface_none_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                         &
         arg_type(GH_FIELD, GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, WTHETA),                   &! gradrinr
         arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! surf_interp
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z0m_eff
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rhostar
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! recip_l_mo_sea
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1_sd
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1) &! q1_sd
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_surface_none_code
  end type bl_surface_none_kernel_type

  public :: bl_surface_none_code

contains

  !> @param[in]     nlayers         Number of layers
  !> @param[in]     wetrho_in_w3    Wet density in w3 space
  !> @param[in,out] rhokm_bl        Momentum eddy diffusivity on BL levels
  !> @param[in,out] rhokh_bl        Heat eddy diffusivity on BL levels
  !> @param[in,out] moist_flux_bl   Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl    Vertical heat flux on BL levels
  !> @param[in,out] gradrinr        Gradient Richardson number in wth
  !> @param[in,out] surf_interp     Surface variables for regridding
  !> @param[in,out] blend_height_tq Blending height for wth levels
  !> @param[in,out] z0m_eff         Grid mean effective roughness length
  !> @param[in,out] ustar           Friction velocity
  !> @param[in,out] rhostar         Surface density
  !> @param[in,out] recip_l_mo_sea  Inverse Obukhov length over sea only
  !> @param[in,out] t1_sd           StDev of level 1 temperature
  !> @param[in,out] q1_sd           StDev of level 1 humidity
  !> @param[in]     ndf_w3          Number of DOFs per cell for w3 space
  !> @param[in]     undf_w3         Number of unique DOFs for w3 space
  !> @param[in]     map_w3          Dofmap for the cell at the base of the column for w3 space
  !> @param[in]     ndf_wth         Number of DOFs per cell for wth space
  !> @param[in]     undf_wth        Number of unique DOFs for wth space
  !> @param[in]     map_wth         Dofmap for the cell at the base of the column for wth space
  !> @param[in]     ndf_surf        Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf       Number of unique DOFs for surface variables
  !> @param[in]     map_surf        Dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_2d          Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d         Number of unique DOFs for 2D fields
  !> @param[in]     map_2d          Dofmap for the cell at the base of the column for 2D fields
  subroutine bl_surface_none_code(nlayers,                       &
                                  wetrho_in_w3,                  &
                                  rhokm_bl,                      &
                                  rhokh_bl,                      &
                                  moist_flux_bl,                 &
                                  heat_flux_bl,                  &
                                  gradrinr,                      &
                                  surf_interp,                   &
                                  blend_height_tq,               &
                                  z0m_eff,                       &
                                  ustar,                         &
                                  rhostar,                       &
                                  recip_l_mo_sea,                &
                                  t1_sd,                         &
                                  q1_sd,                         &
                                  ndf_w3, undf_w3, map_w3,       &
                                  ndf_wth, undf_wth, map_wth,    &
                                  ndf_surf, undf_surf, map_surf, &
                                  ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in), dimension(ndf_surf) :: map_surf
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)   :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in),    dimension(undf_w3)   :: wetrho_in_w3
    real(kind=r_def),    intent(inout), dimension(undf_wth)  :: rhokm_bl
    real(kind=r_def),    intent(inout), dimension(undf_w3)   :: rhokh_bl
    real(kind=r_def),    intent(inout), dimension(undf_w3)   :: moist_flux_bl
    real(kind=r_def),    intent(inout), dimension(undf_w3)   :: heat_flux_bl
    real(kind=r_def),    intent(inout), dimension(undf_wth)  :: gradrinr
    real(kind=r_def),    intent(inout), dimension(undf_surf) :: surf_interp
    integer(kind=i_def), intent(inout), dimension(undf_2d)   :: blend_height_tq
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: z0m_eff
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: ustar
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: rhostar
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: recip_l_mo_sea
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: t1_sd
    real(kind=r_def),    intent(inout), dimension(undf_2d)   :: q1_sd

    ! Surface exchange coefficients and fluxes: free slip, zero flux
    rhokm_bl(map_wth(1))     = 0.0_r_def
    rhokh_bl(map_w3(1))      = 0.0_r_def
    moist_flux_bl(map_w3(1)) = 0.0_r_def
    heat_flux_bl(map_w3(1))  = 0.0_r_def
    gradrinr(map_wth(1))     = 0.0_r_def

    ! 2D variables that need interpolating to cell faces, in the same
    ! slots as set by jules_exp_kernel. Slot 5 (zhnl) is set by bl_exp.
    surf_interp(map_surf(1)+0)  = 0.0_r_def ! flandg: no land
    surf_interp(map_surf(1)+1)  = 0.0_r_def ! rhokm_land
    surf_interp(map_surf(1)+2)  = 0.0_r_def ! rhokm_ssi
    surf_interp(map_surf(1)+3)  = 1.0_r_def ! flandfac
    surf_interp(map_surf(1)+4)  = 1.0_r_def ! fseafac
    surf_interp(map_surf(1)+6)  = 0.0_r_def ! fb_surf
    surf_interp(map_surf(1)+7)  = 0.0_r_def ! k_blend_uv - 1
    surf_interp(map_surf(1)+8)  = 0.0_r_def ! cdr10m
    surf_interp(map_surf(1)+9)  = 0.0_r_def ! cd10m_n
    surf_interp(map_surf(1)+10) = 0.0_r_def ! cdr10m_n

    ! Blending height is used as an array index, so is the lowest level
    blend_height_tq(map_2d(1)) = 1_i_def
    z0m_eff(map_2d(1))         = 0.0_r_def
    ustar(map_2d(1))           = 0.0_r_def

    ! Variables passed to explicit BL
    rhostar(map_2d(1))         = wetrho_in_w3(map_w3(1))
    recip_l_mo_sea(map_2d(1))  = 0.0_r_def
    t1_sd(map_2d(1))           = 0.0_r_def
    q1_sd(map_2d(1))           = 0.0_r_def

  end subroutine bl_surface_none_code

end module bl_surface_none_kernel_mod
