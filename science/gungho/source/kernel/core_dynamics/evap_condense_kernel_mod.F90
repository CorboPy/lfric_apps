!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs a simple warm-rain cloud microphysics scheme with 
!!        latent heating/cooling and basic autoconversion threshold.

!> @details Given the atmospheric temperature and pressure, this kernel computes
!!          the saturation mixing ratio of water vapour. Any excess vapour is
!!          condensed to cloud liquid and split between cloud and rain liquid
!!          according an autoconversion threshold. The rain liquid is then 
!!          advected downwards at a constant speed, with an exponential limiter 
!!          to ensure stability for CFL>1. Afterwards, we evaporate cloud liquid
!!          first before rain liquid, which is reinserted as vapour. The 
!!          potential temperature is adjusted to capture the effects of the 
!!          latent heat release or absorption associated with the phase changes.

module evap_condense_kernel_mod

  use argument_mod,                only: arg_type, GH_SCALAR,                  &
                                         GH_FIELD, GH_WRITE, GH_READ,          &
                                         CELL_COLUMN, GH_REAL
  use constants_mod,               only: r_def, i_def, r_second
  use driver_water_constants_mod,  only: Lv0 => latent_heat_h2o_condensation
  use fs_continuity_mod,           only: Wtheta
  use kernel_mod,                  only: kernel_type
  use physics_common_mod,          only: qsaturation
  use planet_config_mod,           only: recip_epsilon, kappa,                 &
                                         cpd => cp, Rd, p_zero

  implicit none

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: evap_condense_kernel_type
    private
    type(arg_type) :: meta_args(17) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        & ! theta_inc
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! theta_n    
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        & ! mr_v_inc
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        & ! mr_cl_inc
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        & ! mr_r_inc
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_v_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_cl_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_r_n
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! exner_at_wt
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! dz_wtheta
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

  !> @brief Performs a simple warm-rain cloud microphysics scheme
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
  !> @param[in]     cpd          Heat capacity of dry air at const. pressure
  !> @param[in]     cpv          Heat capacity of water vap at const. pressure
  !> @param[in]     cl           Heat capacity of liquid water
  !> @param[in]     p_zero       Reference pressure
  !> @param[in]     dt           The model timestep length
  !> @param[in]     ndf_wtheta   Number of DoF per cell for wtheta
  !> @param[in]     udf_wtheta   Number of total DoF for wtheta
  !> @param[in]     map_wtheta   DoF map for the cell at the base of the column
                 
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
    integer(kind=i_def), intent(in)                           :: nlayers,      & 
                                                                 ndf_wtheta,   & 
                                                                 undf_wtheta
    real(kind=r_def), dimension(undf_wtheta), intent(inout)   :: theta_inc
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: theta_n
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: exner_at_wt
    real(kind=r_def), dimension(undf_wtheta), intent(inout)   :: mr_v_inc,     &
                                                                 mr_cl_inc,    &
                                                                 mr_r_inc      
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: mr_v_n,       &
                                                                 mr_cl_n,      &
                                                                 mr_r_n        
    real(kind=r_def), intent(in)                              :: Rd, Rv, cpd,  &
                                                                 cpv, cl,      &
                                                                 p_zero 
    real(kind=r_second), intent(in)                           :: dt             !> ToDo: figure out fast vs slow dt
    real(kind=r_def), dimension(undf_wtheta), intent(in)      :: dz_wtheta
    integer(kind=i_def), dimension(ndf_wtheta), intent(in)    :: map_wtheta

    ! Internal variables
    real(kind=r_def), dimension(undf_wtheta) :: theta_np1, mr_v_np1,           &
      mr_cl_np1, mr_r_np1, mr_v_np2, mr_cl_np2, mr_r_np2, mr_r_np3
    real(kind=r_def)                         :: mr_sat(0:nlayers),           &
      dm_v_cond(0:nlayers), dm_v_evap(0:nlayers)
    real(kind=r_def)                         :: temperature(0:nlayers),      &
      pressure(0:nlayers)
    real(kind=r_def)                         :: Lv(0:nlayers),               &
      Rm(0:nlayers), cpm(0:nlayers), cvm(0:nlayers)                      
    real(kind=r_def)                         :: rain_carry, rain_available,    &
      fall_fraction, rain_out ! rainout variables
    real(kind=r_def)                         :: cvd, cvv, kappa
    real(kind=r_def), parameter              :: ref_temperature = 273.15_r_def
    real(kind=r_def), parameter              :: cl_threshold = 0.001_r_def      ! cloud threshold
    real(kind=r_def), parameter              :: rain_fall_speed = 5.0_r_def    ! m/s
    integer(kind=i_def)                      :: k

    ! Calculate kappa, cvd, cvv
    kappa = cpd / Rd
    cvd = cpd - Rd
    cvv = cpv - Rv

    ! Convert to temperature and pressure
    do k = 0, nlayers
      temperature(k) = theta_n(map_wtheta(1) + k)                              &
          * exner_at_wt(map_wtheta(1) + k)
      pressure(k) = p_zero * exner_at_wt(map_wtheta(1) + k)                    &
          ** (1.0_r_def/kappa)
    end do

    ! --------------------------------------------------------------------------
    ! Step 1: Calculate saturation mixing ratio
    ! --------------------------------------------------------------------------
    do k = 0, nlayers
      ! This function takes pressure in mbar so divide by 100
      mr_sat(k) = qsaturation(temperature(k), 0.01_r_def*pressure(k))
    end do

    !---------------------------------------------------------------------------
    ! Step 2: Compute change in vapour due to condensation only
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Get latent heat of condensation/evaporation at the current temperature
      Lv(k) = Lv0 - (cl - cpv)*(temperature(k) - ref_temperature)

      ! Compute the change in water vapour mixing ratio
      dm_v_cond(k) = - (mr_v_n(map_wtheta(1) + k) - mr_sat(k)) /               &
              (1.0_r_def + (mr_sat(k) * Lv(k) ** 2.0_r_def) /                  &
                            (cpd * Rv * temperature(k) ** 2.0_r_def))

      ! Clip to prevent evaporation (evap is handled after the rain advection)
      if (dm_v_cond(k) > 0.0_r_def) then 
        dm_v_cond(k) = 0.0_r_def
      end if

      ! Update the water vapour mixing ratio
      mr_v_np1(map_wtheta(1) + k) = mr_v_n(map_wtheta(1) + k) + dm_v_cond(k)
    end do
    ! dm_v_cond is the change in water vapour mixing ratio due to condensation
    ! note that it is negative for condensation (vapour decreases)
  
    !---------------------------------------------------------------------------
    ! Step 3: Partition condensate between cloud and rain
    !       First fill cloud up to a threshold, then split remaining condensate 
    !       between cloud and rain according to cl_threshold.
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Condensation: Add to cloud liquid first
      mr_cl_np1(map_wtheta(1) + k) = mr_cl_n(map_wtheta(1) + k) - dm_v_cond(k)

      ! If cl_threshold is exceeded, immediately convert into rain liquid
      if (mr_cl_np1(map_wtheta(1) + k) > cl_threshold) then
        ! rain update = previous rain + excess cloud
        mr_r_np1(map_wtheta(1) + k)  = mr_r_n(map_wtheta(1) + k)               &
            + (mr_cl_np1(map_wtheta(1) + k) - cl_threshold) 
        ! cloud update = capped at threshold
        mr_cl_np1(map_wtheta(1) + k) = cl_threshold 
      else
        ! rain update = previous rain (no change)
        mr_r_np1(map_wtheta(1) + k)  = mr_r_n(map_wtheta(1) + k) 
      end if
    end do

    !---------------------------------------------------------------------------
    ! Step 4: Advect rain downwards
    !---------------------------------------------------------------------------
    ! Constant-speed rainout using a simple top-down reservoir model.
    ! Rain is passed from the layer above to the layer below with a fixed
    ! fall speed, using an exponential limiter to keep the update stable.
    rain_carry = 0.0_r_def  ! rain falling into current layer from layer above
                            ! starts at 0 for the top layer
    do k = nlayers, 0, -1

      ! add rain carried from layer above to the rain available in this layer
      rain_available = mr_r_np1(map_wtheta(1) + k) + rain_carry   

      ! compute fraction of rain that falls out of this layer in time dt
      fall_fraction = 1.0_r_def - exp(-rain_fall_speed * dt /                  &
          dz_wtheta(map_wtheta(1) + k))

      ! ensure fall_fraction is between 0 and 1
      fall_fraction = max(0.0_r_def, min(fall_fraction, 1.0_r_def))   

      ! compute amount of rain that falls out of this layer
      rain_out = rain_available * fall_fraction
      ! update the rain mixing ratio
      mr_r_np2(map_wtheta(1) + k) = rain_available - rain_out
      rain_carry = rain_out
    end do


    !---------------------------------------------------------------------------
    ! Step 5: Compute change in vapour due to evaporation only
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Compute change in mr_v using *updated* mr_v_np1 and mr_sat from before
      dm_v_evap(k) = - (mr_v_np1(map_wtheta(1) + k) - mr_sat(k)) /             &
              (1.0_r_def + (mr_sat(k) * Lv(k) ** 2.0_r_def) /                  &
                            (cpd * Rv * temperature(k) ** 2.0_r_def))

      ! Clip to prevent condensation (cond is handled before rain advection)
      if (dm_v_evap(k) < 0.0_r_def) then 
        dm_v_evap(k) = 0.0_r_def
      end if

      ! Also clip to prevent evaporation of more than the available condensate
      if (dm_v_evap(k) > (mr_cl_np1(map_wtheta(1) + k)                         & 
            + mr_r_np2(map_wtheta(1) + k))) then
        dm_v_evap(k) = mr_cl_np1(map_wtheta(1) + k)                            &
            + mr_r_np2(map_wtheta(1) + k)
      end if

      ! Update the water vapour mixing ratio
      mr_v_np2(map_wtheta(1) + k) = mr_v_np1(map_wtheta(1) + k) + dm_v_evap(k)
    end do
    ! dm_v_evap is the change in water vapour mixing ratio due to evaporation,
    ! is positive for evaporation (vapour increases) and limited by mr_r+mr_cl

    !---------------------------------------------------------------------------
    ! Step 6: Update cloud and rain liquid after evaporation
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Evaporation: Take from cloud liquid first
      mr_cl_np2(map_wtheta(1) + k) = mr_cl_np1(map_wtheta(1) + k)              & 
            - dm_v_evap(k)  ! mr_cl update = previous cloud - evaporation

      ! Once cloud liquid is depleted, remove from rain liquid
      if (mr_cl_np2(map_wtheta(1) + k) < 0.0_r_def) then
        ! rain update = previous rain + excess cloud (negative)
        mr_r_np3(map_wtheta(1) + k)  = mr_r_np2(map_wtheta(1) + k)             &
            + mr_cl_np2(map_wtheta(1) + k) 
        ! cloud update = capped at zero (it has been fully evaporated)
        mr_cl_np2(map_wtheta(1) + k) = 0.0_r_def 
        ! dm_v_evap cannot evaporate more than the total condensate available 
        ! (due to the check in step 5), so mr_r_np3 *should* not go negative
      else
        ! rain update = previous rain (no change)
        mr_r_np3(map_wtheta(1) + k) = mr_r_np2(map_wtheta(1) + k) 
      end if
    end do

    !---------------------------------------------------------------------------
    ! Step 7: Update the increments for theta, mr_v, mr_cl, and mr_r
    !---------------------------------------------------------------------------
    do k = 0, nlayers
      ! Get the specific heat capacity of the mixture 
      ! (dry air + water vapour + cloud liquid + rain liquid)
      cpm(k) = cpd + mr_v_n(map_wtheta(1) + k) * cpv                           & 
          + (mr_cl_n(map_wtheta(1) + k) + mr_r_n(map_wtheta(1) + k)) * cl      
      cvm(k) = cvd + mr_v_n(map_wtheta(1) + k) * cvv                           &
          + (mr_cl_n(map_wtheta(1) + k) + mr_r_n(map_wtheta(1) + k)) * cl      
      Rm(k) = Rd + mr_v_n(map_wtheta(1) + k) * Rv   
      ! should I be using mr_v_np2, mr_r_np3, mr_cl_np2 here?

      ! Mixing ratios
      mr_v_inc(map_wtheta(1) + k)  = mr_v_np2(map_wtheta(1) + k)               &  
          - mr_v_n(map_wtheta(1) + k)
      mr_cl_inc(map_wtheta(1) + k) = mr_cl_np2(map_wtheta(1) + k)              &  
          - mr_cl_n(map_wtheta(1) + k)
      mr_r_inc(map_wtheta(1) + k)  = mr_r_np3(map_wtheta(1) + k)               &  
          - mr_r_n(map_wtheta(1) + k)

      ! Potential temperature increment due to latent heating
      theta_np1(map_wtheta(1) + k) = theta_n(map_wtheta(1) + k) *              & 
          (1.0_r_def - (dm_v_cond(k) + dm_v_evap(k)) *                         &
          ((cvd * Lv(k) / (cvm(k) * cpd * temperature(k))) -                   &
          (Rv / cvm(k)) * (1.0_r_def - ((Rd * cpm(k)) /                        &
          (cpd * Rm(k))))))

      theta_inc(map_wtheta(1) + k) = theta_np1(map_wtheta(1) + k)              &
          - theta_n(map_wtheta(1) + k)
    end do
    !> ToDo: ensure q_v, q_c, and q_r are above or equal to 0
    !> ToDo: verify total water mass is the same before vs after the scheme

  end subroutine evap_condense_code

end module evap_condense_kernel_mod