!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Collection of routines that are needed for physics.

!> @detail Collection of routines that are needed for physics. These may be
!>         replaced/overloaded by specific schemes as they are brought in from
!>         the UM

module physics_common_mod

  use constants_mod,                 only: r_def
  use planet_config_mod,             only: epsilon
  implicit none
  private

  public qsaturation

contains

  ! Function to return the saturation mr over water
  ! Based on Tetens' formula
  ! es = 6.109 * exp(17.2693882*(T-273.15)/(T-35.86))
  ! qs = epsilon * es / (p - es)
  function qsaturation (T, p) result(qs)
    implicit none
    real(kind=r_def), intent(in) :: T     ! Temperature in Kelvin
    real(kind=r_def), intent(in) :: p     ! Pressure in mb

    real(kind=r_def)             :: qs
    real(kind=r_def)             :: es

    real(kind=r_def),  parameter :: tk0c = 273.15_r_def      ! Temperature of freezing in Kelvin
    ! real(kind=r_def),  parameter :: qsa1 = 3.8_r_def         ! Top constant in qsat equation
    real(kind=r_def),  parameter :: qsa2 = -17.2693882_r_def ! Constant in qsat equation in Kelvin
    real(kind=r_def),  parameter :: qsa3 = 35.86_r_def      ! Constant in qsat equation
    real(kind=r_def),  parameter :: qsa4 = 6.109_r_def       ! Constant in qsat equation in mbar

    if (T > qsa3) then
      es = qsa4 * exp(-qsa2 * (T - tk0c) / (T-qsa3))
      if (p > es) then
        qs = epsilon * es / (p - es)
      else
        qs = 999.0_r_def
      end if
    else
      qs =999.0_r_def
    end if
  end function qsaturation

end module physics_common_mod
