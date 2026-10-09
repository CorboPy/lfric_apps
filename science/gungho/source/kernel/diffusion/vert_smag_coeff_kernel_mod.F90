!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Stability-dependent Smagorinsky eddy diffusivities on theta levels.
!> @details Uses the Lilly (1962) form
!>
!>            K_m = lambda**2 * sqrt( max(0, S**2 - N**2/Pr_t) ),
!>            K_h = K_m / Pr_t,
!>
!>          so mixing stops where N**2/S**2 >= Pr_t and grows in unstable
!>          air. N**2 comes from a parcel test that uses the model's own
!>          thermodynamics: a parcel is moved to a neighbouring level, is
!>          condensed or evaporated there if needed, and its density is
!>          compared with the air already at that level. The test does not
!>          assume dilute vapour, so it treats unsaturated and saturated air,
!>          latent heating, the molecular weight of vapour and condensate
!>          loading in one calculation.
!>
!>          Saturation comes from qsaturation in physics_common_mod, and the
!>          latent heat follows Kirchhoff's law with the same switch
!>          (theta_moist_source) as evap_condense_kernel_mod.
module vert_smag_coeff_kernel_mod

  use argument_mod,               only : arg_type,          &
                                         GH_FIELD, GH_REAL, &
                                         GH_READ, GH_WRITE, &
                                         CELL_COLUMN
  use constants_mod,              only : r_def, i_def
  use driver_water_constants_mod, only : latent_heat_h2o_condensation, &
                                         gas_constant_h2o,             &
                                         heat_capacity_h2o_vapour,     &
                                         heat_capacity_h2o,            &
                                         T_freeze_h2o
  use formulation_config_mod,     only : theta_moist_source
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type
  use mixing_config_mod,          only : smag_prandtl
  use physics_common_mod,         only : qsaturation
  use planet_config_mod,          only : gravity, rd, cp, p_zero, &
                                         one_over_kappa, recip_epsilon

  implicit none

  private

  ! Newton iterations and temperature step for the saturation adjustment
  integer(kind=i_def), parameter :: n_newton = 3
  real(kind=r_def),    parameter :: dt_qsat  = 0.01_r_def

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: vert_smag_coeff_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),  & ! visc_m
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),  & ! visc_h
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),  & ! n_sq
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! theta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! m_t
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! exner
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! shear
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),  & ! lambda_sq
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta)   & ! height_wth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vert_smag_coeff_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vert_smag_coeff_code
  public :: lilly_viscosity

contains

!> @brief Compute N**2 and the eddy diffusivities on theta levels.
!> @param[in]     nlayers    Number of layers in the mesh
!> @param[in,out] visc_m     Momentum eddy diffusivity K_m
!> @param[in,out] visc_h     Scalar eddy diffusivity K_h
!> @param[in,out] n_sq       Squared buoyancy frequency N**2
!> @param[in]     theta      Potential temperature
!> @param[in]     m_v        Vapour mixing ratio
!> @param[in]     m_cl       Cloud liquid mixing ratio
!> @param[in]     m_t        Total water mixing ratio (all species)
!> @param[in]     exner      Exner pressure on theta levels
!> @param[in]     shear      3D strain rate S
!> @param[in]     lambda_sq  Squared mixing length (mix_factor * Delta)**2
!> @param[in]     height_wth Height of theta levels above the surface
!> @param[in]     ndf_wt     Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wt    Number of unique degrees of freedom for Wtheta
!> @param[in]     map_wt     Dofmap for the cell at the base of the column
subroutine vert_smag_coeff_code( nlayers,                &
                                 visc_m, visc_h, n_sq,   &
                                 theta, m_v, m_cl, m_t,  &
                                 exner, shear,           &
                                 lambda_sq, height_wth,  &
                                 ndf_wt, undf_wt, map_wt &
                               )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_wt), intent(inout) :: visc_m, visc_h, n_sq
  real(kind=r_def), dimension(undf_wt), intent(in)    :: theta, m_v, m_cl, m_t
  real(kind=r_def), dimension(undf_wt), intent(in)    :: exner, shear
  real(kind=r_def), dimension(undf_wt), intent(in)    :: lambda_sq, height_wth

  ! Internal variables
  integer(kind=i_def) :: k, b
  real(kind=r_def)    :: temperature(0:nlayers), pressure(0:nlayers)
  real(kind=r_def)    :: t_v(0:nlayers)
  real(kind=r_def)    :: t_v_parcel, n_sq_up, n_sq_dn
  real(kind=r_def)    :: rv_m, cpv_m, cl_m

  ! Moist contributions to the gas constant and heat capacities are only
  ! included with theta_moist_source, as in evap_condense_kernel_mod
  if (theta_moist_source) then
    rv_m  = gas_constant_h2o
    cpv_m = heat_capacity_h2o_vapour
    cl_m  = heat_capacity_h2o
  else
    rv_m  = 0.0_r_def
    cpv_m = 0.0_r_def
    cl_m  = 0.0_r_def
  end if

  b = map_wt(1)

  ! Environment temperature, pressure and virtual temperature. The virtual
  ! temperature is exact for an ideal-gas mixture plus condensate loading,
  ! and matches the moist factors used by the dynamics.
  do k = 0, nlayers
    temperature(k) = theta(b+k) * exner(b+k)
    pressure(k)    = p_zero * exner(b+k)**one_over_kappa
    t_v(k)         = temperature(k) * (1.0_r_def + m_v(b+k)*recip_epsilon) &
                   / (1.0_r_def + m_t(b+k))
  end do

  ! No mixing at the ground or the model top
  visc_m(b) = 0.0_r_def
  visc_h(b) = 0.0_r_def
  n_sq(b)   = 0.0_r_def
  visc_m(b+nlayers) = 0.0_r_def
  visc_h(b+nlayers) = 0.0_r_def
  n_sq(b+nlayers)   = 0.0_r_def

  do k = 1, nlayers - 1

    ! Lift the parcel at level k to level k+1: positive N**2 if it ends up
    ! denser (colder in virtual temperature) than the air there
    t_v_parcel = parcel_virtual_temperature(temperature(k), pressure(k),   &
                                            pressure(k+1), m_v(b+k),       &
                                            m_cl(b+k), m_t(b+k),           &
                                            rv_m, cpv_m, cl_m)
    n_sq_up = gravity * (t_v(k+1) - t_v_parcel)                            &
            / (t_v(k+1) * (height_wth(b+k+1) - height_wth(b+k)))

    ! Lower it to level k-1: positive N**2 if it ends up lighter there.
    ! Level 0 is the ground, which is not part of the mixed column.
    if (k > 1) then
      t_v_parcel = parcel_virtual_temperature(temperature(k), pressure(k), &
                                              pressure(k-1), m_v(b+k),     &
                                              m_cl(b+k), m_t(b+k),         &
                                              rv_m, cpv_m, cl_m)
      n_sq_dn = gravity * (t_v_parcel - t_v(k-1))                          &
              / (t_v(k-1) * (height_wth(b+k) - height_wth(b+k-1)))
      n_sq(b+k) = 0.5_r_def * (n_sq_up + n_sq_dn)
    else
      n_sq(b+k) = n_sq_up
    end if

    visc_m(b+k) = lilly_viscosity(lambda_sq(b+k), shear(b+k), n_sq(b+k), &
                                  smag_prandtl)
    visc_h(b+k) = visc_m(b+k) / smag_prandtl

  end do

end subroutine vert_smag_coeff_code

!> @brief Lilly (1962) momentum eddy diffusivity.
!> @param[in] lambda_sq Squared mixing length
!> @param[in] shear     Strain rate S
!> @param[in] n_sq      Squared buoyancy frequency N**2
!> @param[in] prandtl   Turbulent Prandtl number Pr_t
!> @return    visc_m    lambda**2 * sqrt( max(0, S**2 - N**2/Pr_t) )
pure function lilly_viscosity(lambda_sq, shear, n_sq, prandtl) result(visc_m)

  implicit none

  real(kind=r_def), intent(in) :: lambda_sq, shear, n_sq, prandtl
  real(kind=r_def)             :: visc_m

  visc_m = lambda_sq * sqrt(max(0.0_r_def, shear**2 - n_sq/prandtl))

end function lilly_viscosity

!> @brief Virtual temperature of a parcel moved reversibly from p_a to p_b.
!> @details The parcel first follows the model's unsaturated adiabat. It is
!>          then brought to saturation at constant pressure while conserving
!>          moist enthalpy
!>
!>            H = cp_m(m_v)*T + m_v*l_ref,  l_ref = Lv0 - (cpv - cl)*T0,
!>
!>          which is equivalent to a Kirchhoff latent heat
!>          L(T) = Lv0 + (cpv - cl)*(T - T0). Vapour condenses into cloud if
!>          supersaturated, and cloud (not rain) evaporates if subsaturated.
!>          All condensate stays with the parcel as loading.
!> @param[in] t_a   Parcel temperature at the start
!> @param[in] p_a   Pressure at the start (Pa)
!> @param[in] p_b   Pressure at the end (Pa)
!> @param[in] m_v   Vapour mixing ratio at the start
!> @param[in] m_cl  Cloud liquid mixing ratio at the start
!> @param[in] m_t   Total water mixing ratio (conserved)
!> @param[in] rv_m  Vapour gas constant used in the moist adiabat
!> @param[in] cpv_m Vapour heat capacity used in the moist adiabat
!> @param[in] cl_m  Liquid heat capacity used in the moist adiabat
!> @return    t_v   Parcel virtual temperature at p_b
function parcel_virtual_temperature(t_a, p_a, p_b, m_v, m_cl, m_t, &
                                    rv_m, cpv_m, cl_m) result(t_v)

  implicit none

  real(kind=r_def), intent(in) :: t_a, p_a, p_b, m_v, m_cl, m_t
  real(kind=r_def), intent(in) :: rv_m, cpv_m, cl_m
  real(kind=r_def)             :: t_v

  real(kind=r_def)    :: kappa_m, t_b, p_hpa, q_sat, dqsdt
  real(kind=r_def)    :: l_ref, latent, enthalpy, x
  integer(kind=i_def) :: iter

  ! Unsaturated step along the model's adiabat. kappa_m reduces to rd/cp,
  ! i.e. conserved dry theta, when the moist constants are zero.
  kappa_m = (rd + m_v*rv_m) / cp_moist(m_v)
  t_b     = t_a * (p_b / p_a)**kappa_m

  p_hpa = 0.01_r_def * p_b
  q_sat = qsaturation(t_b, p_hpa)
  x     = m_v

  if ( m_v > q_sat .or. (m_cl > 0.0_r_def .and. m_v < q_sat) ) then

    l_ref    = latent_heat_h2o_condensation - (cpv_m - cl_m) * T_freeze_h2o
    enthalpy = cp_moist(m_v) * t_b + m_v * l_ref

    ! Newton iterations on x - q_sat(T(x)) = 0, where T(x) keeps H fixed so
    ! that dT/dx = -L(T)/cp_m. dq_sat/dT is taken numerically so that any
    ! saturation formula in qsaturation can be used.
    do iter = 1, n_newton
      q_sat  = qsaturation(t_b, p_hpa)
      dqsdt  = ( qsaturation(t_b + dt_qsat, p_hpa)                    &
               - qsaturation(t_b - dt_qsat, p_hpa) ) / (2.0_r_def * dt_qsat)
      latent = l_ref + (cpv_m - cl_m) * t_b
      x      = x - (x - q_sat) / (1.0_r_def + dqsdt * latent / cp_moist(x))
      ! Cannot condense more vapour than there is, or evaporate more cloud
      x      = min(max(x, 0.0_r_def), m_v + m_cl)
      t_b    = (enthalpy - x * l_ref) / cp_moist(x)
    end do

  end if

  t_v = t_b * (1.0_r_def + x*recip_epsilon) / (1.0_r_def + m_t)

contains

  !> Heat capacity of dry air, vapour x and condensate (m_t - x) per kg of
  !> dry air
  pure function cp_moist(x_v) result(cp_m)
    real(kind=r_def), intent(in) :: x_v
    real(kind=r_def)             :: cp_m
    cp_m = cp + x_v*cpv_m + (m_t - x_v)*cl_m
  end function cp_moist

end function parcel_virtual_temperature

end module vert_smag_coeff_kernel_mod
