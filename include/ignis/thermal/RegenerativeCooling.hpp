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
/// `nusselt_multiplier` and `hot_gas_multiplier` make the empirical uncertainty
/// explicit and are dispersed in the Monte Carlo campaign.
///
/// HOT-GAS SIDE
/// ------------
/// The film coefficient comes from the integral boundary layer by default
/// (BoundaryLayer.hpp), or from Bartz on request.  The boundary layer depends
/// on the wall temperature (through its reference state and driving enthalpy)
/// and the wall temperature depends on the film coefficient, so under that
/// model the coupled solve is repeated: march the layer with the current wall
/// temperatures, solve the jacket, update the wall, until the hot-wall
/// temperature stops moving (0.5 K).  Within each pass the local balance still
/// uses the trial wall temperature in the closure, so each station is exactly
/// balanced; only the layer's thickness lags by one pass.
///
/// FILM COOLING
/// ------------
/// A wall film injected at `film_x` (FilmCooling.hpp, Hatch & Papell) lowers
/// the temperature the gas drives the wall towards:
///
///   T_drive = T_aw - eta (T_aw - t_film)
///
/// while h_g is left as the hot-gas model gives it for the unfilmed wall -- the
/// usual superposition, which treats the film as changing the driving
/// temperature and not the layer's conductance.  The film temperature is an
/// input, or by default the jacket's coolant outlet temperature (a fuel film
/// drawn from the injector manifold after the jacket), which couples back
/// through the passes like the boundary layer's wall temperature does.
///
/// TAPERED CHANNELS
/// ----------------
/// `channel_height_profile` and `channel_width_profile` make the channel
/// height and width (fraction of pitch, or absolute in fixed-width mode)
/// tables in x, linear between points and constant beyond them.  The coolant
/// momentum term is then -G du with u = G/rho, which reduces to the
/// constant-area -G^2 d(1/rho) above when the channel does not change.
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

#include <array>
#include <string>
#include <vector>

#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"
#include "ignis/thermal/BoundaryLayer.hpp"
#include "ignis/thermal/CoolantFluid.hpp"
#include "ignis/thermal/FilmCooling.hpp"
#include "ignis/thermal/HeatTransfer.hpp"
#include "ignis/thermo/Transport.hpp"

namespace ignis {

/// Which model supplies the hot-gas film coefficient.
enum class HotGasModel {
  /// Integral turbulent boundary layer marched from the injector face
  /// (BoundaryLayer.hpp).  The default: it carries the layer's history, and it
  /// is closer than Bartz to both measured datasets in validation.md.
  kBoundaryLayer,
  /// Bartz's 1957 closed-form correlation (HeatTransfer.hpp).
  kBartz,
};
std::string toString(HotGasModel m);
HotGasModel hotGasModelFromString(const std::string& s);

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
  HotGasModel hot_gas_model = HotGasModel::kBoundaryLayer;
  /// Explicit scaling of the hot-gas film coefficient, whichever model supplies
  /// it.  The Monte Carlo campaign disperses it; 1 means the model as published.
  double hot_gas_multiplier = 1.0;
  /// Wall temperature assumed upstream of the jacket when the cooled extent
  /// does not start at the injector face (the boundary layer still starts
  /// there).  0 means "the first cooled station's temperature".
  double upstream_wall_temperature = 0.0;
  double gas_emissivity = 0.0;      ///< 0 disables the radiation term
  WallMaterial wall;
  int num_segments = 200;
  double wall_tolerance = 1.0e-6;   ///< K, on the T_wg balance

  /// Channel taper: (x [m], value) tables, linear between points, constant
  /// beyond the ends.  Empty keeps the constant channel_height / width above.
  /// The width table holds the fraction of pitch or, in fixed mode, metres.
  std::vector<std::array<double, 2>> channel_height_profile;
  std::vector<std::array<double, 2>> channel_width_profile;

  /// Wall film (Hatch & Papell).  Off when film_mass_flow is zero.
  double film_mass_flow = 0.0;      ///< kg/s injected along the wall
  /// Slot height, m.  For a ring of film-coolant orifices, the equivalent
  /// slot: total orifice area over the circumference.
  double film_slot_height = 0.5e-3;
  double film_x = 0.0;              ///< m, injection station
  /// Film temperature at injection, K.  0: the jacket's coolant outlet
  /// temperature (fuel film drawn from the manifold after the jacket).
  double film_temperature = 0.0;
  /// Fluid table for the film's properties; empty means the jacket coolant.
  std::string film_coolant;
  double film_pressure = 0.0;       ///< Pa at the slot; 0 = the chamber pressure

  /// Limit on the wetted (coolant-side) wall temperature, K; 0 = no check.
  /// For a hydrocarbon fuel this is its coking limit -- the temperature above
  /// which it deposits carbon on the channel wall.  Ignis ships no default:
  /// see docs/configuration.md for the one figure that could be sourced.
  double coolant_wall_limit = 0.0;
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
  double channel_height = 0.0;   ///< m (varies only with a taper)
  // film (eta = 0 and T_drive = T_aw without one)
  double film_effectiveness = 0.0;
  double t_drive = 0.0;          ///< K, the temperature h_g drives the wall towards
  // boundary layer (zero under the Bartz model)
  double bl_momentum_thickness = 0.0;   ///< m
  double bl_enthalpy_thickness = 0.0;   ///< m
  double bl_re_theta = 0.0;
  int iterations = 0;
  double flux_residual = 0.0; ///< |q_gas - q_cool| / max(|q_gas|, 1 kW/m^2) at convergence
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
  /// The material limit this result was judged against, carried out so that a
  /// reader (or a UI) can draw the constraint without re-reading the library.
  double wall_limit_temperature = 0.0;    ///< K
  std::string wall_material;
  /// |sum q dA - mdot (h_out - h_in)| / |sum q dA|.  Conservation check of the
  /// march itself.
  double energy_balance_residual = 0.0;
  /// Largest per-station |q_gas - q_cool| / q_gas.  Convergence of the coupled
  /// wall balance.
  double max_flux_residual = 0.0;
  /// The model that supplied h_g, and for the boundary layer the number of
  /// wall-temperature passes the coupled solve took to settle.
  std::string hot_gas_model;
  int wall_iterations = 0;
  /// Film, when there is one: what was injected, how fast relative to the gas,
  /// and how far it protects the wall (where eta falls below 0.5 and 0.2; the
  /// jacket end if it does not).
  bool has_film = false;
  double film_mass_flow = 0.0;          ///< kg/s
  double film_temperature = 0.0;        ///< K at injection
  double film_coolant_velocity = 0.0;   ///< m/s
  double film_velocity_ratio = 0.0;     ///< V_g / V_c at the slot
  double film_length_half = 0.0;        ///< m from the slot to eta = 0.5
  double film_length_fifth = 0.0;       ///< m from the slot to eta = 0.2
  bool boiling_detected = false;
  bool wall_limit_exceeded = false;
  /// The hottest wetted wall, and whether it passed `coolant_wall_limit`.
  double max_coolant_side_wall_temperature = 0.0;   ///< K
  double max_coolant_side_wall_temperature_x = 0.0; ///< m
  bool coolant_wall_limit_exceeded = false;
  /// Set only when a *converged* station needed k(T) outside the material's
  /// fitted range.  Root-finder trial evaluations sweep the whole bracket and
  /// are deliberately not recorded here: they would report temperatures the
  /// engine is never predicted to reach.
  bool conductivity_extrapolated = false;
  double conductivity_extrapolation_min = 1.0e30;   ///< K, converged mean wall T
  double conductivity_extrapolation_max = -1.0e30;  ///< K
  double conductivity_extrapolation_x_min = 1.0e30;   ///< m, axial extent
  double conductivity_extrapolation_x_max = -1.0e30;  ///< m
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
