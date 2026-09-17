// SPDX-License-Identifier: MIT
#include "ignis/thermal/RegenerativeCooling.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>

namespace ignis {

std::string toString(ChannelWidthMode m) {
  return m == ChannelWidthMode::kFractionOfPitch ? "fraction_of_pitch" : "fixed";
}

ChannelWidthMode channelWidthModeFromString(const std::string& s) {
  if (s == "fraction_of_pitch") return ChannelWidthMode::kFractionOfPitch;
  if (s == "fixed") return ChannelWidthMode::kFixed;
  throw ConfigError("channel width mode must be 'fraction_of_pitch' or 'fixed', got '" + s + "'");
}

namespace {

/// Colebrook-White friction factor, solved by fixed-point iteration on 1/sqrt(f).
double colebrookFriction(double reynolds, double relative_roughness) {
  if (reynolds < 2300.0) return 64.0 / std::max(reynolds, 1.0);
  double inv_sqrt_f = -2.0 * std::log10(relative_roughness / 3.7 + 5.74 / std::pow(reynolds, 0.9));
  for (int i = 0; i < 100; ++i) {
    const double next = -2.0 * std::log10(relative_roughness / 3.7 +
                                          2.51 * inv_sqrt_f / reynolds);
    if (std::abs(next - inv_sqrt_f) < 1e-12) { inv_sqrt_f = next; break; }
    inv_sqrt_f = next;
  }
  const double f = 1.0 / (inv_sqrt_f * inv_sqrt_f);
  if (reynolds < 4000.0) {
    // Blend laminar and turbulent across the transition band rather than
    // jumping, and flag it to the caller through the low Reynolds number.
    const double w = (reynolds - 2300.0) / 1700.0;
    return (1.0 - w) * (64.0 / reynolds) + w * f;
  }
  return f;
}

/// Nusselt number for the coolant side.
double nusseltNumber(const std::string& correlation, double reynolds, double prandtl,
                     double friction) {
  if (reynolds < 2300.0) return 4.36;  // fully developed laminar, constant q''
  if (correlation == "dittus-boelter") {
    return 0.023 * std::pow(reynolds, 0.8) * std::pow(prandtl, 0.4);
  }
  if (correlation == "gnielinski") {
    const double f8 = friction / 8.0;
    return f8 * (reynolds - 1000.0) * prandtl /
           (1.0 + 12.7 * std::sqrt(f8) * (std::pow(prandtl, 2.0 / 3.0) - 1.0));
  }
  throw ConfigError("cooling: nusselt_correlation must be 'dittus-boelter' or 'gnielinski', got '" +
                    correlation + "'");
}

struct Geometry1D {
  double width = 0.0, land = 0.0, height = 0.0, area = 0.0, perimeter = 0.0, dh = 0.0;
  double pitch = 0.0;
};

Geometry1D channelGeometry(const CoolingSpec& spec, double radius) {
  Geometry1D g;
  g.pitch = 2.0 * constants::pi * radius / spec.num_channels;
  g.height = spec.channel_height;
  if (spec.width_mode == ChannelWidthMode::kFractionOfPitch) {
    g.width = spec.width_fraction * g.pitch;
  } else {
    g.width = spec.channel_width;
  }
  g.land = g.pitch - g.width;
  if (g.land < spec.min_land_width) {
    std::ostringstream os;
    os << "cooling: at r = " << radius * 1e3 << " mm the " << spec.num_channels
       << " channels leave a land of only " << g.land * 1e3 << " mm (minimum "
       << spec.min_land_width * 1e3 << " mm). Reduce the channel count or width.";
    throw ConfigError(os.str());
  }
  g.area = g.width * g.height;
  g.perimeter = 2.0 * (g.width + g.height);
  g.dh = 4.0 * g.area / g.perimeter;
  return g;
}

}  // namespace

std::string CoolingResult::summary() const {
  std::ostringstream os;
  os << std::fixed;
  os << "regenerative cooling\n"
     << "  cooled wetted area     " << std::setprecision(5) << cooled_area << " m^2\n"
     << "  total heat load        " << std::setprecision(2) << total_heat_load * 1e-3 << " kW\n"
     << "  peak heat flux         " << std::setprecision(2) << max_heat_flux * 1e-6
     << " MW/m^2 at x = " << std::setprecision(1) << max_heat_flux_x * 1e3 << " mm\n"
     << "  peak hot-wall temp     " << std::setprecision(1) << max_wall_temperature
     << " K at x = " << max_wall_temperature_x * 1e3 << " mm\n"
     << "  coolant in / out       " << std::setprecision(2) << coolant_inlet_temperature
     << " K -> " << coolant_outlet_temperature << " K  (rise "
     << coolant_temperature_rise << " K)\n"
     << "  coolant pressure drop  " << std::setprecision(4) << coolant_pressure_drop * 1e-6
     << " MPa  (outlet " << coolant_outlet_pressure * 1e-6 << " MPa)\n"
     << "  energy balance residual " << std::scientific << std::setprecision(2)
     << energy_balance_residual << "\n"
     << "  max local flux residual " << max_flux_residual;
  if (!warnings.empty()) {
    os << "\n  warnings:";
    for (const auto& w : warnings) os << "\n    - " << w;
  }
  return os.str();
}

namespace {

/// Shared driver for the cooled and the prescribed-wall surveys.
CoolingResult solveCoolingImpl(const NozzleFlow& flow, const NozzleGeometry& geom,
                               const ChamberResult& chamber, const CoolingSpec& spec,
                               const TransportModel& transport, bool coupled,
                               double prescribed_wall) {
  const double x_lo = spec.x_start;
  const double x_hi = (spec.x_end > 0.0) ? spec.x_end : geom.exitPosition();
  if (!(x_hi > x_lo)) throw ConfigError("cooling: x_end must exceed x_start");
  if (x_lo < 0.0 || x_hi > geom.exitPosition() * (1.0 + 1e-12))
    throw ConfigError("cooling: the cooled extent lies outside the contour");
  if (spec.num_segments < 10) throw ConfigError("cooling: num_segments must be at least 10");
  if (spec.num_channels < 1) throw ConfigError("cooling: num_channels must be at least 1");
  if (!(spec.wall_thickness > 0.0)) throw ConfigError("cooling: wall_thickness must be positive");
  if (coupled) {
    if (!(spec.coolant_mass_flow > 0.0))
      throw ConfigError("cooling: coolant_mass_flow must be positive");
    if (!(spec.inlet_pressure > 0.0))
      throw ConfigError("cooling: inlet_pressure must be positive");
    if (!(spec.inlet_temperature > 0.0))
      throw ConfigError("cooling: inlet_temperature must be positive");
  }

  // --- Bartz reference from the chamber ---------------------------------
  const auto X_chamber = chamber.state.moleFractions();
  const auto tr = transport.mixture(X_chamber, chamber.state.T, chamber.state.cp_frozen);
  BartzReference bref;
  bref.throat_diameter = 2.0 * geom.throatRadius();
  bref.curvature_radius = geom.bartzCurvatureRadius();
  bref.chamber_pressure = chamber.state.p;
  bref.c_star = chamber.c_star;
  bref.chamber_temperature = chamber.state.T;
  bref.viscosity = tr.viscosity;
  bref.cp = chamber.state.cp_frozen;
  bref.prandtl = tr.prandtl;
  bref.multiplier = spec.bartz_multiplier;

  CoolingResult res;
  if (tr.covered_mole_fraction < 0.999) {
    std::ostringstream os;
    os << "transport properties cover only "
       << 100.0 * tr.covered_mole_fraction << " % of the chamber mixture by mole fraction";
    res.warnings.push_back(os.str());
  }

  std::unique_ptr<CoolantFluid> fluid;
  if (coupled) fluid = std::make_unique<CoolantFluid>(CoolantFluid::load(spec.coolant));

  // Segment boundaries, uniform in x.
  const int N = spec.num_segments;
  std::vector<double> edges(static_cast<std::size_t>(N) + 1);
  for (int i = 0; i <= N; ++i) edges[static_cast<std::size_t>(i)] = x_lo + (x_hi - x_lo) * i / N;

  // March direction: counterflow means the coolant enters at x_hi.
  const bool reverse = spec.counterflow;
  const double mdot_total = spec.coolant_mass_flow;
  const double mdot_ch = coupled ? mdot_total / spec.num_channels : 0.0;

  double h_bulk = 0.0, p_bulk = 0.0, T_bulk = 0.0;
  if (coupled) {
    const auto inlet = fluid->at(spec.inlet_temperature, spec.inlet_pressure);
    if (inlet.phase == CoolantPhase::kVapor)
      throw InfeasibleError(
          "cooling: the coolant is already a sub-critical vapour at the jacket inlet (" +
          std::to_string(spec.inlet_temperature) + " K at " +
          std::to_string(spec.inlet_pressure * 1e-6) + " MPa, T_sat = " +
          std::to_string(inlet.t_saturation) + " K); regenerative cooling needs a liquid or "
          "supercritical feed");
    h_bulk = inlet.h;
    p_bulk = spec.inlet_pressure;
    T_bulk = spec.inlet_temperature;
    res.coolant_inlet_temperature = T_bulk;
  }

  std::vector<ThermalStation> stations(static_cast<std::size_t>(N));
  double q_area_sum = 0.0;
  double rho_prev = 0.0;
  bool was_liquid = false;
  bool near_sat_warned = false;

  for (int k = 0; k < N; ++k) {
    const int i = reverse ? (N - 1 - k) : k;
    const double xa = edges[static_cast<std::size_t>(i)];
    const double xb = edges[static_cast<std::size_t>(i) + 1];
    const double xm = 0.5 * (xa + xb);
    const double r = geom.radius(xm);
    const double ra = geom.radius(xa), rb = geom.radius(xb);
    const double dA_gas = constants::pi * (ra + rb) * std::hypot(xb - xa, rb - ra);
    const double dl = std::hypot(xb - xa, rb - ra);   // wetted length along the wall

    ExpansionState gas;
    const double eps = std::max(geom.area(xm) / geom.throatArea(), 1.0);
    gas = flow.atAreaRatio(eps, xm > geom.throatPosition());

    ThermalStation st;
    st.x = xm;
    st.radius = r;
    st.area_ratio = eps;
    st.gas_T = gas.gas.T;
    st.gas_p = gas.gas.p;
    st.gas_mach = gas.mach;
    st.gas_gamma = gas.gas.gamma_s;
    st.t_adiabatic_wall = recoveryTemperature(gas.gas.T, gas.mach, gas.gas.gamma_s, bref.prandtl);

    Geometry1D cg{};
    double A_cool_per_len = 0.0;
    CoolantState cs{};
    if (coupled) {
      cg = channelGeometry(spec, r);
      st.channel_width = cg.width;
      st.land_width = cg.land;
      st.hydraulic_diameter = cg.dh;
      cs = fluid->at(T_bulk, p_bulk);
      // Boiling is only meaningful if the coolant was a sub-critical liquid and
      // then crossed its saturation line.  A supercritical fluid, or one above
      // its critical temperature, has no liquid branch to boil out of.
      if (cs.phase == CoolantPhase::kLiquid) was_liquid = true;
      if (was_liquid && cs.phase == CoolantPhase::kVapor && !res.boiling_detected) {
        res.boiling_detected = true;
        std::ostringstream os;
        os << "coolant boiling at x = " << xm * 1e3 << " mm: the bulk crossed T_sat = "
           << cs.t_saturation << " K at " << p_bulk * 1e-6
           << " MPa; the single-phase correlations no longer apply";
        res.warnings.push_back(os.str());
      }
      if (cs.near_saturation && !near_sat_warned) {
        near_sat_warned = true;
        std::ostringstream os;
        os << "coolant within " << 100.0 * CoolantState::kSaturationMargin
           << " % of its saturation temperature at x = " << xm * 1e3
           << " mm (T = " << T_bulk << " K, T_sat = " << cs.t_saturation << " K)";
        res.warnings.push_back(os.str());
      }
      st.coolant_phase = toString(cs.phase);
      st.coolant_T = T_bulk;
      st.coolant_p = p_bulk;
      st.coolant_rho = cs.rho;
      st.coolant_cp = cs.cp;
      st.coolant_mu = cs.mu;
      st.coolant_k = cs.k;
      st.coolant_h = cs.h;
      st.mass_flux = mdot_ch / cg.area;
      st.coolant_velocity = st.mass_flux / cs.rho;
      st.reynolds = st.mass_flux * cg.dh / cs.mu;
      st.prandtl = cs.prandtl;
      st.friction_factor = colebrookFriction(st.reynolds, spec.roughness / cg.dh);
      st.nusselt = spec.nusselt_multiplier *
                   nusseltNumber(spec.nusselt_correlation, st.reynolds, st.prandtl,
                                 st.friction_factor);
      st.h_coolant = st.nusselt * cs.k / cg.dh;
    }

    // --- solve the coupled wall balance ---------------------------------
    auto wallConductivity = [&](double t_hot, double t_cold) {
      const double tm = 0.5 * (t_hot + t_cold);
      if (!spec.wall.inValidRange(tm)) res.conductivity_extrapolated = true;
      return spec.wall.conductivity(tm);
    };

    // residual(T_wg) = q_gas - q_cool, both per m^2 of gas-side wall.
    auto residual = [&](double t_wg, ThermalStation* out) {
      const double hg = bartzFilmCoefficient(bref, eps, gas.mach, gas.gas.gamma_s, t_wg);
      const double qc = hg * (st.t_adiabatic_wall - t_wg);
      const double qr = grayGasRadiation(gas.gas.T, t_wg, spec.gas_emissivity,
                                         spec.wall.emissivity);
      const double q = qc + qr;
      // Through-wall conduction: iterate once on the mean-temperature k.
      double t_wc = t_wg;
      double kw = spec.wall.conductivity(t_wg);
      for (int it = 0; it < 30; ++it) {
        const double R = cylindricalWallResistance(r, spec.wall_thickness, kw);
        const double t_new = t_wg - q * R;
        const double kn = wallConductivity(t_wg, t_new);
        if (std::abs(t_new - t_wc) < 1e-10 && std::abs(kn - kw) < 1e-10) { t_wc = t_new; kw = kn; break; }
        t_wc = t_new;
        kw = kn;
      }
      if (out != nullptr) {
        out->h_gas = hg;
        out->q_convective = qc;
        out->q_radiative = qr;
        out->q_total = q;
        out->t_wall_hot = t_wg;
        out->t_wall_cold = t_wc;
        out->wall_conductivity = kw;
      }
      if (!coupled) return 0.0;
      // Fin-augmented coolant-side conductance referred to the gas-side area.
      const double m = std::sqrt(2.0 * st.h_coolant / (kw * cg.land));
      const double mh = m * cg.height;
      const double eta = (mh > 1e-12) ? std::tanh(mh) / mh : 1.0;
      A_cool_per_len = cg.width + 2.0 * cg.height * eta;
      if (out != nullptr) out->fin_efficiency = eta;
      const double ratio = A_cool_per_len / cg.pitch;   // per unit gas-side area
      const double q_cool = st.h_coolant * ratio * (t_wc - T_bulk);
      return q - q_cool;
    };

    if (!coupled) {
      residual(prescribed_wall, &st);
      st.iterations = 1;
      st.flux_residual = 0.0;
    } else {
      // Bracket: T_wg between the coolant bulk temperature and T_aw.
      double lo = T_bulk + 1.0e-3;
      double hi = st.t_adiabatic_wall - 1.0e-6;
      double f_lo = residual(lo, nullptr);
      double f_hi = residual(hi, nullptr);
      if (f_lo * f_hi > 0.0) {
        std::ostringstream os;
        os << "cooling: the wall energy balance is not bracketed at x = " << xm * 1e3
           << " mm (residuals " << f_lo << " and " << f_hi << " W/m^2 between " << lo
           << " K and " << hi << " K)";
        throw ConvergenceError(os.str());
      }
      int it = 0;
      double t_wg = 0.5 * (lo + hi);
      for (; it < 300; ++it) {
        t_wg = 0.5 * (lo + hi);
        const double f = residual(t_wg, nullptr);
        if (f * f_lo > 0.0) { lo = t_wg; f_lo = f; } else { hi = t_wg; f_hi = f; }
        if (hi - lo < spec.wall_tolerance) break;
      }
      const double f_final = residual(t_wg, &st);
      st.iterations = it + 1;
      st.flux_residual = (st.q_total != 0.0) ? std::abs(f_final / st.q_total) : 0.0;
      res.max_flux_residual = std::max(res.max_flux_residual, st.flux_residual);
    }

    if (st.t_wall_hot > spec.wall.max_temperature) res.wall_limit_exceeded = true;

    st.segment_heat = st.q_total * dA_gas;
    q_area_sum += st.segment_heat;
    res.cooled_area += dA_gas;
    if (st.q_total > res.max_heat_flux) { res.max_heat_flux = st.q_total; res.max_heat_flux_x = xm; }
    if (st.t_wall_hot > res.max_wall_temperature) {
      res.max_wall_temperature = st.t_wall_hot;
      res.max_wall_temperature_x = xm;
    }

    if (coupled) {
      // Energy: march the bulk enthalpy so the balance closes identically.
      h_bulk += st.segment_heat / mdot_total;
      // Momentum + friction pressure change.
      const double G = st.mass_flux;
      const double dp_fric = st.friction_factor * (dl / cg.dh) * (G * G / (2.0 * cs.rho));
      double dp_mom = 0.0;
      if (rho_prev > 0.0) dp_mom = G * G * (1.0 / cs.rho - 1.0 / rho_prev);
      rho_prev = cs.rho;
      p_bulk -= (dp_fric + dp_mom);
      if (!(p_bulk > 0.0)) {
        std::ostringstream os;
        os << "cooling: the coolant pressure reached zero at x = " << xm * 1e3
           << " mm; the channel cannot pass " << mdot_total << " kg/s";
        throw InfeasibleError(os.str());
      }
      T_bulk = fluid->temperatureFromEnthalpy(h_bulk, p_bulk, T_bulk);
    }
    stations[static_cast<std::size_t>(i)] = st;
  }

  res.stations = std::move(stations);
  res.total_heat_load = q_area_sum;
  if (coupled) {
    res.coolant_outlet_temperature = T_bulk;
    res.coolant_temperature_rise = T_bulk - res.coolant_inlet_temperature;
    res.coolant_outlet_pressure = p_bulk;
    res.coolant_pressure_drop = spec.inlet_pressure - p_bulk;
    const auto in = fluid->at(res.coolant_inlet_temperature, spec.inlet_pressure);
    const auto out = fluid->at(T_bulk, p_bulk);
    const double q_coolant = mdot_total * (out.h - in.h);
    res.energy_balance_residual =
        (q_area_sum != 0.0) ? std::abs(q_area_sum - q_coolant) / std::abs(q_area_sum) : 0.0;
    if (res.wall_limit_exceeded) {
      std::ostringstream os;
      os << "peak hot-wall temperature " << res.max_wall_temperature
         << " K exceeds the " << spec.wall.max_temperature << " K limit of "
         << spec.wall.name;
      res.warnings.push_back(os.str());
    }
    if (res.conductivity_extrapolated)
      res.warnings.push_back("wall conductivity was evaluated outside its fitted temperature range");
  }
  return res;
}

}  // namespace

CoolingResult solveRegenerativeCooling(const NozzleFlow& flow, const NozzleGeometry& geom,
                                       const ChamberResult& chamber, const CoolingSpec& spec,
                                       const TransportModel& transport) {
  return solveCoolingImpl(flow, geom, chamber, spec, transport, true, 0.0);
}

CoolingResult surveyHotGasSide(const NozzleFlow& flow, const NozzleGeometry& geom,
                               const ChamberResult& chamber, const CoolingSpec& spec,
                               const TransportModel& transport,
                               double prescribed_wall_temperature) {
  if (!(prescribed_wall_temperature > 0.0))
    throw ConfigError("hot-gas survey: the prescribed wall temperature must be positive");
  return solveCoolingImpl(flow, geom, chamber, spec, transport, false,
                          prescribed_wall_temperature);
}

}  // namespace ignis
