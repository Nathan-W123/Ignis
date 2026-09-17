// SPDX-License-Identifier: MIT
#pragma once
/// \file RegenerativeCooling.hpp
/// \brief One-dimensional regenerative cooling-channel model.
///
/// MODEL
/// -----
/// The cooled length is divided into axial segments.  In each segment the
/// hot gas, the wall and the coolant are coupled through a single unknown, the
/// hot-side wall temperature T_wg, by requiring that the same heat crosses all
/// three resistances:
///
///   q_gas  = h_g (T_aw - T_wg) + q_rad(T_wg)                [per m^2 of gas-side wall]
///   T_wc   = T_wg - q_gas R_wall(T_mean)                    [radial conduction]
///   q_cool = h_c (A_cool_eff / A_gas) (T_wc - T_bulk)       [fin-augmented convection]
///
/// and solving q_gas(T_wg) = q_cool(T_wg) to tolerance.  Because h_g depends on
/// T_wg through Bartz's sigma factor and h_c depends on T_wc through the
/// coolant properties, the balance is genuinely nonlinear and is solved with a
/// bracketed secant/bisection hybrid.
///
/// The land between two channels is treated as a straight rectangular fin of
/// height equal to the channel height:
///
///   m = sqrt(2 h_c / (k_wall t_land)),  eta_fin = tanh(m H) / (m H)
///   A_cool_eff / length = w_channel + 2 H eta_fin
///
/// The outer closeout is taken as adiabatic, which is conservative.
///
/// The coolant is marched in *enthalpy*, not temperature, so the energy
/// balance closes identically regardless of how violently cp varies near the
/// pseudo-critical line.  Pressure is marched with both the friction term and
/// the momentum (acceleration) term, which is significant when a supercritical
/// coolant expands by a factor of five along the jacket:
///
///   dp = -f (dx / D_h) (G^2 / (2 rho)) - G^2 d(1/rho)
///
/// CORRELATIONS (all empirical, all labelled)
/// ------------------------------------------
///   Dittus-Boelter   Nu = 0.023 Re^0.8 Pr^0.4
///     F. W. Dittus and L. M. K. Boelter, Univ. Calif. Publ. Eng. 2, 443 (1930).
///   Gnielinski       Nu = (f/8)(Re-1000) Pr / [1 + 12.7 sqrt(f/8)(Pr^(2/3)-1)]
///     V. Gnielinski, Int. Chem. Eng. 16, 359-368 (1976).
///   Colebrook-White  1/sqrt(f) = -2 log10[eps/(3.7 D) + 2.51/(Re sqrt f)]
///     C. F. Colebrook, J. Inst. Civ. Eng. 11, 133-156 (1939).
///   Laminar          f = 64/Re,  Nu = 4.36 (constant heat flux, circular)
///
/// `nusselt_multiplier` and `bartz_multiplier` make the empirical uncertainty
/// explicit and are dispersed in the Monte Carlo campaign.
///
/// FAILURE MODES REPORTED
/// ----------------------
///   * coolant boiling (sub-critical pressure with T_bulk >= T_sat)
///   * coolant state outside the tabulated property range
///   * wall temperature above the material limit
///   * conductivity evaluated outside its fitted range
///   * channel geometry that leaves no land between channels
///   * coolant pressure driven to zero
/// None of these are silently absorbed: each sets a flag, records a message,
/// and (for the hard ones) throws.

#include <string>
#include <vector>

#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"
#include "ignis/thermal/CoolantFluid.hpp"
#include "ignis/thermal/HeatTransfer.hpp"
#include "ignis/thermo/Transport.hpp"

namespace ignis {

/// How the channel width is derived from the local circumference.
enum class ChannelWidthMode {
  kFractionOfPitch,  ///< w = width_fraction * (2 pi r / N); land takes the rest
  kFixed,            ///< w is constant; the land width then varies with r
};
std::string toString(ChannelWidthMode m);
ChannelWidthMode channelWidthModeFromString(const std::string& s);

/// Cooling-jacket definition.
struct CoolingSpec {
  std::string coolant = "methane";
  int num_channels = 100;
  ChannelWidthMode width_mode = ChannelWidthMode::kFractionOfPitch;
  double width_fraction = 0.5;      ///< channel width / channel pitch
  double channel_width = 0.0;       ///< m, used when width_mode == kFixed
  double channel_height = 3.0e-3;   ///< m
  double min_land_width = 0.5e-3;   ///< m, geometric validity floor
  double wall_thickness = 1.0e-3;   ///< m
  double roughness = 5.0e-6;        ///< m, absolute channel roughness
  double coolant_mass_flow = 0.0;   ///< kg/s through the whole jacket
  double inlet_temperature = 0.0;   ///< K
  double inlet_pressure = 0.0;      ///< Pa
  bool counterflow = true;          ///< coolant enters at the nozzle exit
  double x_start = 0.0;             ///< cooled extent, m (injector side)
  double x_end = -1.0;              ///< cooled extent, m (-1 => nozzle exit)
  /// Alternative to `x_end`: end the jacket where the divergent reaches this
  /// area ratio.  Real engines regeneratively cool the chamber, throat and the
  /// first part of the divergent and then switch to a radiation-cooled or
  /// film-cooled extension, and expressing the extent as an area ratio keeps
  /// that split meaningful while the expansion ratio is being swept or
  /// optimised.  Values <= 1 disable it.  The uncooled extension is NOT
  /// modelled: its wall temperature is not predicted and its (small) heat load
  /// is not counted.
  double x_end_area_ratio = 0.0;
  std::string nusselt_correlation = "dittus-boelter";
  double nusselt_multiplier = 1.0;
  double bartz_multiplier = 1.0;
  double gas_emissivity = 0.0;      ///< 0 disables the radiation term
  WallMaterial wall;
  int num_segments = 200;
  double wall_tolerance = 1.0e-6;   ///< K, on the T_wg balance
};

/// Converged state of one axial segment.
struct ThermalStation {
  double x = 0.0;             ///< m
  double radius = 0.0;        ///< m
  double area_ratio = 0.0;
  // gas side
  double gas_T = 0.0, gas_p = 0.0, gas_mach = 0.0, gas_gamma = 0.0;
  double t_adiabatic_wall = 0.0;
  double h_gas = 0.0;         ///< W/(m^2 K)
  double q_convective = 0.0;  ///< W/m^2
  double q_radiative = 0.0;   ///< W/m^2
  double q_total = 0.0;       ///< W/m^2
  // wall
  double t_wall_hot = 0.0, t_wall_cold = 0.0;
  double wall_conductivity = 0.0;
  // coolant
  double coolant_T = 0.0, coolant_p = 0.0, coolant_rho = 0.0, coolant_cp = 0.0;
  double coolant_mu = 0.0, coolant_k = 0.0, coolant_h = 0.0;
  double coolant_velocity = 0.0, mass_flux = 0.0;
  std::string coolant_phase;
  double reynolds = 0.0, prandtl = 0.0, nusselt = 0.0, h_coolant = 0.0;
  double friction_factor = 0.0, fin_efficiency = 0.0;
  double channel_width = 0.0, land_width = 0.0, hydraulic_diameter = 0.0;
  double segment_heat = 0.0;  ///< W absorbed by the coolant in this segment
  int iterations = 0;
  double flux_residual = 0.0; ///< |q_gas - q_cool| / q_gas at convergence
};

/// Result of a cooling analysis.
struct CoolingResult {
  std::vector<ThermalStation> stations;   ///< ordered by increasing x
  double max_wall_temperature = 0.0;      ///< K
  double max_wall_temperature_x = 0.0;    ///< m
  double max_heat_flux = 0.0;             ///< W/m^2
  double max_heat_flux_x = 0.0;           ///< m
  double total_heat_load = 0.0;           ///< W
  double coolant_inlet_temperature = 0.0;
  double coolant_outlet_temperature = 0.0;
  double coolant_temperature_rise = 0.0;  ///< K
  double coolant_pressure_drop = 0.0;     ///< Pa
  double coolant_outlet_pressure = 0.0;   ///< Pa
  double cooled_area = 0.0;               ///< m^2
  /// |sum q dA - mdot (h_out - h_in)| / |sum q dA|.  Conservation check of the
  /// march itself.
  double energy_balance_residual = 0.0;
  /// Largest per-station |q_gas - q_cool| / q_gas.  Convergence of the coupled
  /// wall balance.
  double max_flux_residual = 0.0;
  bool boiling_detected = false;
  bool wall_limit_exceeded = false;
  bool conductivity_extrapolated = false;
  double conductivity_extrapolation_min = 1.0e30;
  double conductivity_extrapolation_max = -1.0e30;
  std::vector<std::string> warnings;
  std::string summary() const;
};

/// Solve the coupled hot-gas / wall / coolant problem.
///
/// \param flow    converged expansion for the operating point
/// \param geom    nozzle contour
/// \param chamber chamber solution (supplies c*, T_c and transport reference)
/// \param spec    jacket definition
/// \param transport transport model over the active species set
CoolingResult solveRegenerativeCooling(const NozzleFlow& flow, const NozzleGeometry& geom,
                                       const ChamberResult& chamber, const CoolingSpec& spec,
                                       const TransportModel& transport);

/// Hot-gas-side-only thermal survey with a prescribed wall temperature.
/// Useful for uncooled or radiation-cooled sections and for verification.
CoolingResult surveyHotGasSide(const NozzleFlow& flow, const NozzleGeometry& geom,
                               const ChamberResult& chamber, const CoolingSpec& spec,
                               const TransportModel& transport, double prescribed_wall_temperature);

}  // namespace ignis
