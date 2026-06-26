!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs a simple condensation/evaporation scheme with latent heating
!!        and basic Rotunno & Emanuel (1987) rainout/re-evaporation.

!> @details Given the atmospheric temperature and pressure, this kernel computes
!!          the saturation mixing ratio of water vapour. Any excess vapour is
!!          condensed to cloud liquid, while any cloud liquid in an unsaturated
!!          environment is evaporated to water vapour. The potential temperature
!!          is adjusted to capture the effects of the latent heat release or
!!          absorption associated with this phase change.
!!          Note: this only works with the lowest order spaces

module evap_condense_kernel_mod

  use argument_mod,                only: arg_type, GH_SCALAR,         &
                                         GH_FIELD, GH_WRITE, GH_READ, &
                                         DOF, GH_REAL
  use constants_mod,               only: r_def, i_def, r_second
  use driver_water_constants_mod,  only: Lv0 => latent_heat_h2o_condensation
  use fs_continuity_mod,           only: Wtheta
  use kernel_mod,                  only: kernel_type
  use physics_common_mod,          only: qsaturation
  use planet_config_mod,           only: recip_epsilon, kappa, cpd => cp, Rd, p_zero

  implicit none

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: evap_condense_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        & ! theta_inc
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! theta_n    
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        & ! mr_v_inc
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        & ! mr_cl_inc
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        & ! mr_r_inc
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! mr_v_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! mr_cl_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! mr_r_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! exner_at_wt
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        & ! dz_wtheta      !> ToDo: pass through fast_physics_alg
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! Rd
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! Rv
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cpd
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cpv
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cl
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! p_zero
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  & ! dt
    /)
    integer :: operates_on = CELL_COLUMN

  contains
      procedure, nopass :: evap_condense_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: evap_condense_code

contains

  !> @brief Performs a simple condensation/evaporation scheme with latent heating
  !> @param[in]     nlayers      Integer the number of layers
  !> @param[in,out] theta_inc    Potential temperature increment
  !> @param[in]     theta_n      Potential temperature input
  !> @param[in,out] mr_v_inc     Water vapour mixing ratio increment
  !> @param[in,out] mr_cl_inc    Liquid cloud mixing ratio increment
  !> @param[in,out] mr_r_inc     Rain liquid mixing ratio increment
  !> @param[in]     mr_v_n       Water vapour mixing ratio input
  !> @param[in]     mr_cl_n      Liquid cloud mixing ratio input
  !> @param[in]     mr_r_n       Rain liquid mixing ratio input
  !> @param[in]     exner_at_wt  Exner pressure at Wtheta points
  !> @param[in]     dz_wtheta    Vertical spacing at Wtheta points
  !> @param[in]     Rd           Gas constant for dry air
  !> @param[in]     Rv           Gas constant for water vapour
  !> @param[in]     cpd          Heat capacity of dry air at constant pressure
  !> @param[in]     cpv          Heat capacity of water vap at constant pressure
  !> @param[in]     cl           Heat capacity of liquid water
  !> @param[in]     p_zero       Reference pressure
  !> @param[in]     dt           The model timestep length
  !> @param[in]     ndf_wtheta    The number of degrees of freedom per cell for wtheta
  !> @param[in]     udf_wtheta    The number of total degrees of freedom for wtheta
  !> @param[in]     map_wtheta    Integer array holding the dofmap for the cell at the base of the column
                 
  subroutine evap_condense_code( nlayers,                                      &
                                 theta_inc,                                    &
                                 theta_n,                                      &
                                 mr_v_inc,                                     &
                                 mr_cl_inc,                                    &
                                 mr_r_inc,                                     &
                                 mr_v_n,                                       & 
                                 mr_cl_n,                                      &
                                 mr_r_n,                                       &
                                 exner_at_wt,                                  &
                                 dz_wtheta,                                    &
                                 Rd,                                           &
                                 Rv,                                           &
                                 cpd,                                          &
                                 cpv,                                          &
                                 cl,                                           &
                                 p_zero,                                       &
                                 dt,                                           &
                                 ndf_wtheta, undf_wtheta, map_wtheta )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                           :: nlayers, ndf_wtheta, undf_wtheta
    real(kind=r_def), dimension(undf_wtheta), intent(inout)   :: theta_inc
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: theta_n
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: exner_at_wt
    real(kind=r_def), dimension(undf_wtheta), intent(inout)   :: mr_v_inc, mr_cl_inc, mr_r_inc
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: mr_v_n, mr_cl_n, mr_r_n
    real(kind=r_def), intent(in)                              :: Rd, Rv, cpd, cpv, cl, p_zero
    real(kind=r_second), intent(in)                           :: dt
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: dz_wtheta
    integer(kind=i_def), dimension(ndf_wtheta), intent(in)    :: map_wtheta

    ! Internal variables
    real(kind=r_def), dimension(undf_wtheta) :: theta_np1, mr_v_np1, mr_cl_np1, mr_r_np1
    real(kind=r_def),                        :: mr_sat(0:nlayers-1), dm_v(0:nlayers-1)
    real(kind=r_def)                         :: temperature(0:nlayers-1), pressure(0:nlayers-1)
    real(kind=r_def),                        :: Lv(0:nlayers-1), Rm(0:nlayers-1), cpm(0:nlayers-1), cvm(0:nlayers-1)
    real(kind=r_def) :: cvd, cvv, kappa
    real(kind=r_def), parameter :: ref_temperature = 273.15_r_def
    real(kind=r_def), parameter :: cl_threshold = 0.001_r_def      ! cloud threshold
    real(kind=i_def) :: k

    ! Convert to temperature and pressure
    kappa = cpd / Rd
    do k = 0, nlayers
      temperature(k) = theta_n(map_wtheta(1) + k) * exner_at_wt(map_wtheta(1) + k)
      pressure(k) = p_zero * exner_at_wt(map_wtheta(1) + k) ** (1.0_r_def/kappa)
    end do
    ! Internal variables don't need map_wtheta as they are only used within the column and not written back to the global field

    ! Thermodynamic quantities
    kappa = Rd / cpd ! float
    cvd = cpd - Rd   ! float
    cvv = cpv - Rv   ! float
    do k = 0, nlayers
      cpm(k) = cpd + mr_v_n(map_wtheta(1) + k) * cpv + (mr_cl_n(map_wtheta(1) + k) + mr_r_n(map_wtheta(1) + k)) * cl     ! should I move this to theta adjustment? (after the updates to mr_v, mr_cl, mr_r)
      cvm(k) = cvd + mr_v_n(map_wtheta(1) + k) * cvv + (mr_cl_n(map_wtheta(1) + k) + mr_r_n(map_wtheta(1) + k)) * cl     ! should I move this to theta adjustment? (after the updates to mr_v, mr_cl, mr_r)
      Rm(k) = Rd + mr_v_n(map_wtheta(1) + k) * Rv                                                                        ! should I move this to theta adjustment? (after the updates to mr_v, mr_cl, mr_r)
      Lv(k) = Lv0 - (cl - cpv)*(temperature(k) - ref_temperature)                                                        ! This is used to calculate dm_v
    end do
    !!> ToDo: break up the above into multiple lines  

    ! Calculate saturation mixing ratio for column
    do k = 0, nlayers
      ! This function takes pressure in mbar so divide by 100
      mr_sat(k) = qsaturation(temperature(k), 0.01_r_def*pressure(k))
    end do

    !---------------------------------------------------------------------------
    ! Compute net vapour change dm_v
    !   dm_v < 0 : condensation (vapour -> condensate) — we lose vapour
    !   dm_v > 0 : evaporation (condensate -> vapour) — we gain vapour
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      dm_v(k) = - (mr_v_n(map_wtheta(1) + k) - mr_sat(k)) /                      &
              (1.0_r_def + (mr_sat(k) * Lv(k) ** 2.0_r_def) /                  &
                            (cpd * Rv * temperature(k) ** 2.0_r_def))
    end do
    ! dm_v represents the change in water vapour mixing ratio due to condensation/evaporation for the column

    ! Clip to prevent evaporating more condensate than available 
    do k = 0, nlayers
      if (dm_v(k) > 0.0_r_def) then                                                       ! if evaporation is happening 
        dm_v(k) = min(dm_v(k), mr_cl_n(map_wtheta(1) + k) + mr_r_n(map_wtheta(1) + k))    ! Limit dm_v gain to the total condensate available (cloud + rain)
      end if
    end do

    ! Provisional vapour after phase change
    do k = 0, nlayers
      mr_v_np1(map_wtheta(1) + k) = mr_v_n(map_wtheta(1) + k) + dm_v(k)
    end do

    !---------------------------------------------------------------------------
    ! Rotunno & Emanuel–style partitioning of condensate
    !
    !   - If dm_v < 0: condensation (vapour -> condensate) — we gain condensate. 
    !       First fill cloud up to a threshold, then split remaining condensate 
    !       between cloud and rain according to cl_threshold.
    !
    !   - If dm_v > 0: evaporation (condensate -> vapour) — we lose condensate. 
    !       Evaporate first from cloud, then from rain.
    !       This is a simple approach to account for the fact that raindrops are
    !       much larger than cloud droplets, and thus more likely to evaporate.
    !       A more comprehensive timescale-based approach could be implemented 
    !       in the future.
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Condensation: Add to cloud liquid first
      ! Evaporation: Remove from cloud first, (then from rain, because raindrops are much larger than cloud droplets)
      mr_cl_np1(map_wtheta(1) + k) = mr_cl_n(map_wtheta(1) + k) - dm_v(k)  ! cloud liquid cell update = previous cloud +/- condensation/evaporation. If dm_v < 0, cloud increases.
      !mr_r_np1(map_wtheta(1) + k)  = mr_r_n(map_wtheta(1) + k)            ! rain liquid cell update = previous rain (updated later)

      if (dm_v(k) < 0.0_r_def) then     ! Condensation

        ! If cloud liquid exceeds threshold, it is immediately converted into rain liquid
        if (mr_cl_np1(map_wtheta(1) + k) > cl_threshold) then
          mr_r_np1(map_wtheta(1) + k)  = mr_r_n(map_wtheta(1) + k) + (mr_cl_np1(map_wtheta(1) + k) - cl_threshold) ! rain update = previous rain + excess cloud
          mr_cl_np1(map_wtheta(1) + k) = cl_threshold ! cloud update = capped at threshold
        end if
      
      else                              ! Evaporation

        ! Once cloud liquid is depleted, remove from rain liquid
        if (mr_cl_np1(map_wtheta(1) + k) < 0.0_r_def) then
          mr_r_np1(map_wtheta(1) + k)  = mr_r_n(map_wtheta(1) + k) + mr_cl_np1(map_wtheta(1) + k) ! rain update = previous rain + excess cloud (negative)
          mr_cl_np1(map_wtheta(1) + k) = 0.0_r_def ! cloud update = capped at zero (it has been fully evaporated)
          ! dm_v cannot evaporate more than the total condensate available, so mr_r_np1 *should* not go negative
        end if
      end if
    end do
    ! for the above: do I need a check for dm_v(k) = 0.0_r_def? (no condensation/evaporation, so no change to cloud/rain)

    ! REST OF SCHEME HERE (theta adjustment, ...then downward advection of rain?)


    ! The rest below this is from the original saturation adjustment scheme

    ! Update fields
    mr_v_np1 = mr_v_n - dm_v
    mr_cl_np1 = mr_cl_n + dm_v
    theta_np1 = theta_n * (                                                    &
        1.0_r_def + dm_v * ((cvd * Lv / (cvm * cpd * temperature))             &
        - (Rv / cvm) * (1 - ((Rd * cpm) / (cpd * Rm))))                        &
    )

    ! Compute final increments !> ToDo PROBABLY INCORRECT AT THE MOMENT! - dm_v is wrong. It should be + dm_v
    theta_inc = theta_np1 - theta_n
    mr_v_inc = - dm_v
    mr_cl_inc = mr_cl_np1 - mr_cl_n

  end subroutine evap_condense_code

end module evap_condense_kernel_mod
