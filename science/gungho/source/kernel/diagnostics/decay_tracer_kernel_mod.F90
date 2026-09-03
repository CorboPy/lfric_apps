!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies the fixed surface concentration to the decaying tracer
!>
!> @details Overwrites the levels of the decaying tracer field nearest the
!>          surface with a prescribed concentration. This provides the fixed
!>          lower boundary condition which acts as the source of the tracer,
!>          the tracer elsewhere being removed by exponential decay.
!>          The tracer is held in Wtheta (full levels), for which level 0 is
!>          the ground, so a reset_level of 1 fixes the concentration at the
!>          surface itself.
module decay_tracer_kernel_mod

  use argument_mod,      only : arg_type,                &
                                GH_FIELD, GH_REAL,       &
                                GH_READ, GH_READWRITE,   &
                                CELL_COLUMN, GH_SCALAR,  &
                                GH_INTEGER
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: decay_tracer_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                         &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta), &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),              &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)               &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: decay_tracer_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: decay_tracer_code

contains

!> @brief Applies the fixed surface concentration to the decaying tracer
!> @param[in]     nlayers       Number of layers
!> @param[in,out] decay_tracer  Decaying tracer field
!> @param[in]     surface_value Concentration to hold the tracer at near the
!>                              surface
!> @param[in]     reset_level   Number of full levels, counting up from the
!>                              ground, held at surface_value
!> @param[in]     ndf_wt        Number of degrees of freedom per cell
!> @param[in]     undf_wt       Total number of degrees of freedom
!> @param[in]     map_wt        Dofmap for the cell at the base of the column
subroutine decay_tracer_code( nlayers,                              &
                              decay_tracer,                         &
                              surface_value,                        &
                              reset_level,                          &
                              ndf_wt,                               &
                              undf_wt,                              &
                              map_wt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt
  integer(kind=i_def), intent(in) :: undf_wt
  real(kind=r_def), dimension(undf_wt), intent(inout) :: decay_tracer
  real(kind=r_def), intent(in) :: surface_value
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in) :: reset_level

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: df

  ! Wtheta holds nlayers+1 levels, indexed 0 (the ground) to nlayers
  do k = 0, min(reset_level, nlayers+1)-1
    do df = 1, ndf_wt
      decay_tracer( map_wt(df) + k ) = surface_value
    end do
  end do

end subroutine decay_tracer_code

end module decay_tracer_kernel_mod
