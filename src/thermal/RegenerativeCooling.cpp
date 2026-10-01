// SPDX-License-Identifier: MIT
#include "ignis/thermal/RegenerativeCooling.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
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

std::string toString(HotGasModel m) {
  return m == HotGasModel::kBartz ? "bartz" : "boundary_layer";
}

HotGasModel hotGasModelFromString(const std::string& s) {
  if (s == "boundary_layer") return HotGasModel::kBoundaryLayer;
  if (s == "bartz") return HotGasModel::kBartz;
  throw ConfigError("hot_gas_model must be 'boundary_layer' or 'bartz', got '" + s + "'");
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

/// Linear interpolation in an (x, value) table, constant beyond its ends.
double profileAt(const std::vector<std::array<double, 2>>& table, double x) {
  if (x <= table.front()[0]) return table.front()[1];
  if (x >= table.back()[0]) return table.back()[1];
  for (std::size_t i = 1; i < table.size(); ++i) {
    if (x <= table[i][0]) {
      const double t = (x - table[i - 1][0]) / (table[i][0] - table[i - 1][0]);
      return table[i - 1][1] + t * (table[i][1] - table[i - 1][1]);
    }
  }
  return table.back()[1];
}

void checkProfile(const std::vector<std::array<double, 2>>& table, const char* name) {
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (!(table[i][1] > 0.0))
      throw ConfigError(std::string("cooling: every value in ") + name + " must be positive");
    if (i > 0 && !(table[i][0] > table[i - 1][0]))
      throw ConfigError(std::string("cooling: the x values of ") + name +
                        " must be strictly increasing");
  }
}

Geometry1D channelGeometry(const CoolingSpec& spec, double radius, double x) {
  Geometry1D g;
  g.pitch = 2.0 * constants::pi * radius / spec.num_channels;
  g.height = spec.channel_height_profile.empty() ? spec.channel_height
                                                 : profileAt(spec.channel_height_profile, x);
  const bool tapered_width = !spec.channel_width_profile.empty();
  if (spec.width_mode == ChannelWidthMode::kFractionOfPitch) {
    const double f = tapered_width ? profileAt(spec.channel_width_profile, x) : spec.width_fraction;
    if (!(f > 0.0 && f < 1.0))
      throw ConfigError("cooling: the channel width fraction must lie between 0 and 1");
    g.width = f * g.pitch;
  } else {
    g.width = tapered_width ? profileAt(spec.channel_width_profile, x) : spec.channel_width;
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
     << "  hot-gas model          " << hot_gas_model;
  if (wall_iterations > 1) os << "  (" << wall_iterations << " coupled passes)";
  os << "\n";
  if (has_film) {
    os << "  wall film              " << std::setprecision(3) << film_mass_flow << " kg/s at "
       << std::setprecision(1) << film_temperature << " K, V_g/V_c " << std::setprecision(2)
       << film_velocity_ratio << ";  eta > 0.5 for " << std::setprecision(1)
       << film_length_half * 1e3 << " mm, > 0.2 for " << film_length_fifth * 1e3 << " mm\n";
  }
  os << "  energy balance residual " << std::scientific << std::setprecision(2)
     << energy_balance_residual << "\n"
     << "  max local flux residual " << max_flux_residual;
  if (!warnings.empty()) {
    os << "\n  warnings:";
    for (const auto& w : warnings) os << "\n    - " << w;
  }
  return os.str();
}

namespace {

/// Gas-side state at one jacket segment, computed once and shared by every
/// pass of the coupled solve.
struct SegmentGas {
  double xa = 0.0, xb = 0.0, xm = 0.0, r = 0.0, dA = 0.0, dl = 0.0, eps = 0.0;
  double t_aw = 0.0;
  ExpansionState gas;
};

/// Shared driver for the cooled and the prescribed-wall surveys.
CoolingResult solveCoolingImpl(const NozzleFlow& flow, const NozzleGeometry& geom,
                               const ChamberResult& chamber, const CoolingSpec& spec,
                               const TransportModel& transport, bool coupled,
                               double prescribed_wall) {
  const double x_lo = spec.x_start;
  double x_hi = (spec.x_end > 0.0) ? spec.x_end : geom.exitPosition();
  if (spec.x_end_area_ratio > 1.0) {
    if (spec.x_end_area_ratio > geom.expansionRatio()) {
      x_hi = geom.exitPosition();
    } else {
      // Walk the divergent to the first station at or beyond the requested area
      // ratio; the contour is monotone there so this is unambiguous.
      x_hi = geom.exitPosition();
      for (const auto& st : geom.stations()) {
        if (st.x <= geom.throatPosition()) continue;
        if (st.area_ratio >= spec.x_end_area_ratio) { x_hi = st.x; break; }
      }
    }
  }
  if (!(x_hi > x_lo)) throw ConfigError("cooling: x_end must exceed x_start");
  if (x_lo < 0.0 || x_hi > geom.exitPosition() * (1.0 + 1e-12))
    throw ConfigError("cooling: the cooled extent lies outside the contour");
  if (spec.x_end > 0.0 && spec.x_end_area_ratio > 1.0)
    throw ConfigError("cooling: give either x_end or x_end_area_ratio, not both");
  if (spec.num_segments < 10) throw ConfigError("cooling: num_segments must be at least 10");
  if (spec.num_channels < 1) throw ConfigError("cooling: num_channels must be at least 1");
  if (!(spec.wall_thickness > 0.0)) throw ConfigError("cooling: wall_thickness must be positive");
  if (!(spec.hot_gas_multiplier > 0.0))
    throw ConfigError("cooling: hot_gas_multiplier must be positive");
  checkProfile(spec.channel_height_profile, "channel_height_profile");
  checkProfile(spec.channel_width_profile, "channel_width_profile");
  if (spec.film_mass_flow < 0.0) throw ConfigError("cooling: film_mass_flow cannot be negative");
  if (spec.film_mass_flow > 0.0 && !(spec.film_slot_height > 0.0))
    throw ConfigError("cooling: film_slot_height must be positive");
  if (coupled) {
    if (!(spec.coolant_mass_flow > 0.0))
      throw ConfigError("cooling: coolant_mass_flow must be positive");
    if (!(spec.inlet_pressure > 0.0))
      throw ConfigError("cooling: inlet_pressure must be positive");
    if (!(spec.inlet_temperature > 0.0))
      throw ConfigError("cooling: inlet_temperature must be positive");
  }

  // --- Bartz reference from the chamber ---------------------------------
  // Also the recovery-factor Prandtl number for both models, so the driving
  // temperature is defined the same way whichever supplies h_g.
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
  bref.multiplier = spec.hot_gas_multiplier;

  CoolingResult base;
  base.hot_gas_model = toString(spec.hot_gas_model);
  if (tr.covered_mole_fraction < 0.999) {
    std::ostringstream os;
    os << "transport properties cover only "
       << 100.0 * tr.covered_mole_fraction << " % of the chamber mixture by mole fraction";
    base.warnings.push_back(os.str());
  }

  std::unique_ptr<CoolantFluid> fluid;
  if (coupled) fluid = std::make_unique<CoolantFluid>(CoolantFluid::load(spec.coolant));

  // Segment boundaries.  Uniform in x, except where the jacket starts at the
  // injector face: the boundary layer starts there, its heat-transfer
  // coefficient varies on the scale of its virtual-origin length (a few mm),
  // and a uniform grid would resolve that peak only at first order in the
  // segment length.  There the first segment is 10/N of the uniform width and
  // each one grows by 1.2 until it reaches it, so the first segment shrinks
  // as 1/N^2.  The grid is the same whichever hot-gas model runs, so the two
  // can be compared station by station.
  const int N = spec.num_segments;
  std::vector<double> edges(static_cast<std::size_t>(N) + 1);
  {
    std::vector<double> width(static_cast<std::size_t>(N), 1.0);
    if (x_lo <= 1.0e-9 * geom.exitPosition()) {
      double w = 1.0 / N;
      for (std::size_t k = 0; k < width.size() && w < 1.0; ++k, w *= 1.5) width[k] = w;
    }
    double total = 0.0;
    for (double w : width) total += w;
    edges[0] = x_lo;
    for (int i = 0; i < N; ++i)
      edges[static_cast<std::size_t>(i) + 1] =
          edges[static_cast<std::size_t>(i)] + (x_hi - x_lo) * width[static_cast<std::size_t>(i)] / total;
    edges[static_cast<std::size_t>(N)] = x_hi;
  }

  // --- gas side at every segment, once -----------------------------------
  auto gasAt = [&](double xa, double xb) {
    SegmentGas s;
    s.xa = xa;
    s.xb = xb;
    s.xm = 0.5 * (xa + xb);
    s.r = geom.radius(s.xm);
    const double ra = geom.radius(xa), rb = geom.radius(xb);
    s.dA = constants::pi * (ra + rb) * std::hypot(xb - xa, rb - ra);
    s.dl = std::hypot(xb - xa, rb - ra);   // wetted length along the wall
    s.eps = std::max(geom.area(s.xm) / geom.throatArea(), 1.0);
    s.gas = flow.atAreaRatio(s.eps, s.xm > geom.throatPosition());
    s.t_aw = recoveryTemperature(s.gas.gas.T, s.gas.mach, s.gas.gas.gamma_s, bref.prandtl);
    return s;
  };
  std::vector<SegmentGas> seg(static_cast<std::size_t>(N));
  for (int i = 0; i < N; ++i)
    seg[static_cast<std::size_t>(i)] =
        gasAt(edges[static_cast<std::size_t>(i)], edges[static_cast<std::size_t>(i) + 1]);

  // --- wall film ----------------------------------------------------------
  // drive[i] is the temperature h_g drives the wall towards at segment i: the
  // recovery temperature, or below it where a film protects the wall.
  std::vector<double> drive(static_cast<std::size_t>(N)), film_eta(static_cast<std::size_t>(N), 0.0);
  for (int i = 0; i < N; ++i) drive[static_cast<std::size_t>(i)] = seg[static_cast<std::size_t>(i)].t_aw;
  const bool film = spec.film_mass_flow > 0.0;
  const bool film_auto_temperature = film && !(spec.film_temperature > 0.0);
  if (film_auto_temperature && !coupled)
    throw ConfigError("film cooling: a gas-side survey has no coolant outlet to draw the film "
                      "from; give film_temperature");
  std::unique_ptr<CoolantFluid> film_fluid;
  std::vector<ExpansionState> film_grid_gas;
  std::vector<double> film_grid_slope;
  struct FilmInfo {
    double coolant_velocity = 0.0, velocity_ratio = 0.0, length_half = 0.0, length_fifth = 0.0;
    bool clamped = false, reaches_half = false, reaches_fifth = false;
  } film_info;
  // Film state for a given injection temperature.  H, the heat-capacity group,
  // is integrated on its own fine grid from the slot, so the slot can sit
  // anywhere (upstream of the jacket too) and the segment grid does not enter.
  auto updateFilm = [&](double t_c) {
    if (!film_fluid)
      film_fluid = std::make_unique<CoolantFluid>(
          CoolantFluid::load(spec.film_coolant.empty() ? spec.coolant : spec.film_coolant));
    const double xf = spec.film_x;
    if (xf < 0.0 || !(xf < x_hi))
      throw ConfigError("film cooling: film_x must lie between the injector face and the jacket end");
    auto gasAtX = [&](double x) {
      const double eps = std::max(geom.area(x) / geom.throatArea(), 1.0);
      return flow.atAreaRatio(eps, x > geom.throatPosition());
    };
    const ExpansionState slot_gas = gasAtX(xf);
    const double p_slot = spec.film_pressure > 0.0 ? spec.film_pressure : slot_gas.gas.p;
    const auto cs = film_fluid->at(t_c, p_slot);
    FilmSlot slot;
    slot.mass_flow = spec.film_mass_flow;
    slot.slot_height = spec.film_slot_height;
    slot.perimeter = 2.0 * constants::pi * geom.radius(xf);
    slot.temperature = t_c;
    slot.density = cs.rho;
    slot.cp = cs.cp;
    slot.conductivity = cs.k;
    slot.gas_velocity = slot_gas.u;
    const double factor = hatchPapellFactor(slot, &film_info.clamped);
    film_info.coolant_velocity = slot.coolantVelocity();
    film_info.velocity_ratio = slot.gas_velocity / film_info.coolant_velocity;

    // h of the report's assumption 6 at (t_g + t_c)/2, in the local core gas.
    // The fine grid's gas states and wall slopes do not depend on the film
    // temperature, so they are computed on the first call only.
    const GasMixture& gm = flow.mixture();
    constexpr int kFine = 300;
    if (film_grid_gas.empty()) {
      for (int k = 0; k <= kFine; ++k) {
        const double x = xf + (x_hi - xf) * k / kFine;
        film_grid_gas.push_back(gasAtX(x));
        const double dr = 1.0e-6 * geom.exitPosition();
        film_grid_slope.push_back((geom.radius(std::min(x + dr, geom.exitPosition())) -
                                   geom.radius(std::max(x - dr, 0.0))) / (2.0 * dr));
      }
    }
    auto integrand = [&](int k, double x) {
      const ExpansionState& g = film_grid_gas[static_cast<std::size_t>(k)];
      const double tf = 0.5 * (g.gas.T + t_c);
      const double cp = gm.cpFrozen(g.gas.n, tf);
      const auto tf_props = transport.mixture(g.gas.moleFractions(), tf, cp);
      const double rho = g.gas.p * g.gas.M / (constants::R_universal * tf);
      const double r = geom.radius(x);
      const double h = hatchPapellGasCoefficient(rho, std::max(g.u, 1.0e-3), 2.0 * r,
                                                 tf_props.viscosity, tf_props.conductivity,
                                                 tf_props.prandtl);
      const double slope = film_grid_slope[static_cast<std::size_t>(k)];
      return h * 2.0 * constants::pi * r * std::sqrt(1.0 + slope * slope) /
             (spec.film_mass_flow * cs.cp);
    };
    std::vector<double> xs(kFine + 1), H(kFine + 1, 0.0);
    double f_prev = 0.0;
    for (int k = 0; k <= kFine; ++k) {
      xs[static_cast<std::size_t>(k)] = xf + (x_hi - xf) * k / kFine;
      const double f_now = integrand(k, xs[static_cast<std::size_t>(k)]);
      if (k > 0)
        H[static_cast<std::size_t>(k)] = H[static_cast<std::size_t>(k) - 1] +
            0.5 * (f_prev + f_now) * (xs[static_cast<std::size_t>(k)] - xs[static_cast<std::size_t>(k) - 1]);
      f_prev = f_now;
    }
    auto etaAt = [&](double x) {
      if (x <= xf) return 0.0;   // no film upstream of the slot
      const double t = (x - xf) / (x_hi - xf) * kFine;
      const std::size_t k = std::min(static_cast<std::size_t>(t), static_cast<std::size_t>(kFine - 1));
      const double w = t - static_cast<double>(k);
      return hatchPapellEffectiveness(H[k] + w * (H[k + 1] - H[k]), factor);
    };
    film_info.reaches_half = film_info.reaches_fifth = false;
    film_info.length_half = film_info.length_fifth = x_hi - xf;
    for (int k = 1; k <= kFine; ++k) {
      const double e = hatchPapellEffectiveness(H[static_cast<std::size_t>(k)], factor);
      const double d = xs[static_cast<std::size_t>(k)] - xf;
      if (!film_info.reaches_half && e < 0.5) { film_info.reaches_half = true; film_info.length_half = d; }
      if (!film_info.reaches_fifth && e < 0.2) { film_info.reaches_fifth = true; film_info.length_fifth = d; }
    }
    for (int i = 0; i < N; ++i) {
      const auto& sg = seg[static_cast<std::size_t>(i)];
      const double e = etaAt(sg.xm);
      film_eta[static_cast<std::size_t>(i)] = e;
      drive[static_cast<std::size_t>(i)] = sg.t_aw - e * (sg.t_aw - t_c);
    }
  };

  // March direction: counterflow means the coolant enters at x_hi.
  const bool reverse = spec.counterflow;
  const double mdot_total = spec.coolant_mass_flow;
  const double mdot_ch = coupled ? mdot_total / spec.num_channels : 0.0;

  // One pass of the coolant march with a given film-coefficient model.
  // hg(i, T_wg) is the hot-gas coefficient at segment i for a trial wall
  // temperature.  Each pass starts from `base`, so warnings are not repeated.
  using FilmModel = std::function<double(int, double)>;
  auto runPass = [&](const FilmModel& hg) {
    CoolingResult res = base;
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
    double rho_prev = 0.0, G_prev = 0.0;
    bool was_liquid = false;
    bool near_sat_warned = false;

    for (int k = 0; k < N; ++k) {
      const int i = reverse ? (N - 1 - k) : k;
      const SegmentGas& sg = seg[static_cast<std::size_t>(i)];
      const double xm = sg.xm;
      const double r = sg.r;
      const ExpansionState& gas = sg.gas;
      const double eps = sg.eps;

      ThermalStation st;
      st.x = xm;
      st.radius = r;
      st.area_ratio = eps;
      st.gas_T = gas.gas.T;
      st.gas_p = gas.gas.p;
      st.gas_mach = gas.mach;
      st.gas_gamma = gas.gas.gamma_s;
      st.t_adiabatic_wall = sg.t_aw;
      st.t_drive = drive[static_cast<std::size_t>(i)];
      st.film_effectiveness = film_eta[static_cast<std::size_t>(i)];

      Geometry1D cg{};
      double A_cool_per_len = 0.0;
      CoolantState cs{};
      if (coupled) {
        cg = channelGeometry(spec, r, xm);
        st.channel_width = cg.width;
        st.channel_height = cg.height;
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
      // The root finder evaluates the wall balance across its whole bracket --
      // from the coolant bulk temperature up to the adiabatic wall temperature --
      // so most evaluations are trial states that are nowhere near the answer.
      // Only the converged evaluation (the one that fills a ThermalStation) is
      // allowed to record an extrapolation, because only that one describes a
      // state the engine is actually predicted to be in.
      // The through-wall fixed point below has trial iterates of its own, so
      // the record is made once, from the state it settles on.
      auto recordConductivity = [&](double t_hot, double t_cold) {
        const double tm = 0.5 * (t_hot + t_cold);
        if (spec.wall.inValidRange(tm)) return;
        res.conductivity_extrapolated = true;
        res.conductivity_extrapolation_min = std::min(res.conductivity_extrapolation_min, tm);
        res.conductivity_extrapolation_max = std::max(res.conductivity_extrapolation_max, tm);
        res.conductivity_extrapolation_x_min = std::min(res.conductivity_extrapolation_x_min, xm);
        res.conductivity_extrapolation_x_max = std::max(res.conductivity_extrapolation_x_max, xm);
      };
      auto wallConductivity = [&](double t_hot, double t_cold) {
        return spec.wall.conductivity(0.5 * (t_hot + t_cold));
      };

      // residual(T_wg) = q_gas - q_cool, both per m^2 of gas-side wall.
      auto residual = [&](double t_wg, ThermalStation* out) {
        const double h_g = hg(i, t_wg);
        const double qc = h_g * (st.t_drive - t_wg);
        const double qr = grayGasRadiation(gas.gas.T, t_wg, spec.gas_emissivity,
                                           spec.wall.emissivity);
        const double q = qc + qr;
        // Through-wall conduction: iterate once on the mean-temperature k.
        // The conductivity depends on the mean wall temperature, which depends on
        // the flux, which depends on the conductivity.  The map is a strong
        // contraction (k varies by a few percent across the wall), so a handful of
        // fixed-point steps to 1e-8 relative is far more than enough.
        double t_wc = t_wg;
        double kw = spec.wall.conductivity(t_wg);
        for (int it = 0; it < 20; ++it) {
          const double R = cylindricalWallResistance(r, spec.wall_thickness, kw);
          const double t_new = t_wg - q * R;
          const double kn = wallConductivity(t_wg, t_new);
          const bool done = std::abs(t_new - t_wc) < 1e-8 * std::max(1.0, t_new);
          t_wc = t_new;
          kw = kn;
          if (done) break;
        }
        if (out != nullptr) {
          recordConductivity(t_wg, t_wc);
          out->h_gas = h_g;
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
        // Bracket.  The residual falls strictly with T_wg (the gas gives less
        // heat to a hotter wall, the coolant takes more), so any interval
        // reaching below both the coolant and the driving temperature and up
        // to the recovery temperature contains exactly one root.  A film
        // colder than the coolant reverses those two, and a film drawn from
        // the jacket outlet can sit within millikelvin of the coolant beside
        // its slot, so the bracket is opened a kelvin beyond them rather than
        // pinched between them.
        double lo = std::min(T_bulk, st.t_drive) - 1.0;
        double hi = std::max(T_bulk + 1.0, st.t_adiabatic_wall - 1.0e-6);
        double f_lo = residual(lo, nullptr);
        double f_hi = residual(hi, nullptr);
        if (f_lo * f_hi > 0.0) {
          std::ostringstream os;
          os << "cooling: the wall energy balance is not bracketed at x = " << xm * 1e3
             << " mm (residuals " << f_lo << " and " << f_hi << " W/m^2 between " << lo
             << " K and " << hi << " K)";
          throw ConvergenceError(os.str());
        }
        // Illinois regula falsi: keeps the bracket but converges super-linearly,
        // so the coupled balance needs roughly ten flux evaluations instead of the
        // thirty a pure bisection to 1e-6 K would take.
        int it = 0;
        double a = lo, b = hi, fa = f_lo, fb = f_hi;
        double t_wg = 0.5 * (a + b);
        for (; it < 200; ++it) {
          if (fb == fa) break;
          t_wg = b - fb * (b - a) / (fb - fa);
          const double guard = 0.01 * (b - a);
          if (!(t_wg > a + guard && t_wg < b - guard)) t_wg = 0.5 * (a + b);
          const double f = residual(t_wg, nullptr);
          if (std::abs(f) < 1e-9 * std::max(1.0, std::abs(st.t_adiabatic_wall) * 1.0e3) ||
              (b - a) < spec.wall_tolerance)
            break;
          if (f * fb < 0.0) { a = b; fa = fb; } else { fa *= 0.5; }
          b = t_wg;
          fb = f;
          if (a > b) { std::swap(a, b); std::swap(fa, fb); }
        }
        const double f_final = residual(t_wg, &st);
        st.iterations = it + 1;
        // Relative to the local flux, floored at 1 kW/m^2: beside a film
        // drawn from the jacket the flux can be nearly zero, and a ratio to
        // zero says nothing about convergence.
        st.flux_residual = std::abs(f_final) / std::max(std::abs(st.q_total), 1.0e3);
        res.max_flux_residual = std::max(res.max_flux_residual, st.flux_residual);
      }

      res.wall_limit_temperature = spec.wall.max_temperature;
      res.wall_material = spec.wall.name;
      if (st.t_wall_hot > spec.wall.max_temperature) res.wall_limit_exceeded = true;

      st.segment_heat = st.q_total * sg.dA;
      q_area_sum += st.segment_heat;
      res.cooled_area += sg.dA;
      if (st.q_total > res.max_heat_flux) { res.max_heat_flux = st.q_total; res.max_heat_flux_x = xm; }
      if (st.t_wall_hot > res.max_wall_temperature) {
        res.max_wall_temperature = st.t_wall_hot;
        res.max_wall_temperature_x = xm;
      }
      if (coupled && st.t_wall_cold > res.max_coolant_side_wall_temperature) {
        res.max_coolant_side_wall_temperature = st.t_wall_cold;
        res.max_coolant_side_wall_temperature_x = xm;
      }

      if (coupled) {
        // Energy: march the bulk enthalpy so the balance closes identically.
        h_bulk += st.segment_heat / mdot_total;
        // Momentum + friction pressure change.  dp = -G du with u = G/rho, so
        // a channel whose area changes accelerates its coolant too; with a
        // constant area this is the usual G^2 (1/rho - 1/rho_prev).
        const double G = st.mass_flux;
        const double dp_fric = st.friction_factor * (sg.dl / cg.dh) * (G * G / (2.0 * cs.rho));
        double dp_mom = 0.0;
        if (rho_prev > 0.0) dp_mom = 0.5 * (G + G_prev) * (G / cs.rho - G_prev / rho_prev);
        rho_prev = cs.rho;
        G_prev = G;
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
      if (spec.coolant_wall_limit > 0.0 &&
          res.max_coolant_side_wall_temperature > spec.coolant_wall_limit) {
        res.coolant_wall_limit_exceeded = true;
        std::ostringstream os;
        os << "coolant-side wall reaches " << res.max_coolant_side_wall_temperature
           << " K at x = " << res.max_coolant_side_wall_temperature_x * 1e3 << " mm, above the "
           << spec.coolant_wall_limit << " K limit set for " << spec.coolant
           << " (for a hydrocarbon fuel, its coking limit)";
        res.warnings.push_back(os.str());
      }
      if (res.conductivity_extrapolated) {
        // Say which end of the fit was left, and where: a cold coolant inlet
        // below the fit is routine, a hot wall above it is not, and one
        // message must not describe the second as the first.
        struct Side { int count = 0; double t_lo = 1e30, t_hi = -1e30, x_lo = 1e30, x_hi = -1e30; };
        Side cold, hot;
        bool peak_flux_inside = true;
        for (const auto& s : res.stations) {
          const double tm = 0.5 * (s.t_wall_hot + s.t_wall_cold);
          if (spec.wall.inValidRange(tm)) continue;
          Side& side = tm < spec.wall.valid_min ? cold : hot;
          ++side.count;
          side.t_lo = std::min(side.t_lo, tm);
          side.t_hi = std::max(side.t_hi, tm);
          side.x_lo = std::min(side.x_lo, s.x);
          side.x_hi = std::max(side.x_hi, s.x);
          if (s.x == res.max_heat_flux_x) peak_flux_inside = false;
        }
        std::ostringstream os;
        os << "wall conductivity was extrapolated outside the fitted range ["
           << spec.wall.valid_min << ", " << spec.wall.valid_max << "] K for "
           << spec.wall.name << " (k(T) is clamped to the fit endpoints there):";
        auto describe = [&](const Side& side, const char* label) {
          os << " " << label << ", converged mean wall temperature " << side.t_lo << " to "
             << side.t_hi << " K over x = " << side.x_lo * 1e3 << " to " << side.x_hi * 1e3
             << " mm (" << side.count << " stations);";
        };
        if (cold.count > 0) describe(cold, "below the fit");
        if (hot.count > 0) describe(hot, "ABOVE the fit");
        if (hot.count == 0)
          os << " only the cold end is affected, which is the inlet of a cryogenic jacket";
        else
          os << " the wall is hotter than the conductivity data extend";
        os << (peak_flux_inside ? "; the peak-flux station is inside the fitted range."
                                : "; the peak-flux station is OUTSIDE the fitted range.");
        res.warnings.push_back(os.str());
      }
    }
    return res;
  };

  const bool use_bl = spec.hot_gas_model == HotGasModel::kBoundaryLayer;

  // --- integral boundary layer -------------------------------------------
  // The layer grows from (effectively) zero thickness at the injector face
  // whatever the cooled extent.  Its origin is a station of its own at x = 0:
  // the jacket's first station is a segment midpoint, and starting the march
  // there would hand that segment the arbitrary starting thickness.  The
  // stretch upstream of the jacket (if any) is added as extra stations.
  std::vector<SegmentGas> pad;
  if (use_bl) {
    pad.push_back(gasAt(0.0, 0.0));
    if (x_lo > 1.0e-9 * geom.exitPosition()) {
      const int P = std::max(8, static_cast<int>(std::ceil(N * x_lo / (x_hi - x_lo))));
      for (int k = 0; k < P; ++k) pad.push_back(gasAt(x_lo * k / P, x_lo * (k + 1) / P));
    }
  }
  const std::size_t n_pad = pad.size();
  const std::size_t n_edge = use_bl ? n_pad + static_cast<std::size_t>(N) : 0;
  auto edgeGas = [&](std::size_t j) -> const SegmentGas& {
    return j < n_pad ? pad[j] : seg[j - n_pad];
  };
  const GasMixture& mix = flow.mixture();
  auto referenceAt = [&](std::size_t j, double T) {
    const ExpansionState& e = edgeGas(j).gas;
    ReferenceState rs;
    rs.rho = e.gas.p * e.gas.M / (constants::R_universal * T);
    rs.cp = mix.cpFrozen(e.gas.n, T);
    const auto t = transport.mixture(e.gas.moleFractions(), T, rs.cp);
    rs.viscosity = t.viscosity;
    rs.prandtl = t.prandtl;
    return rs;
  };
  std::vector<BoundaryLayerEdge> edge(n_edge);
  for (std::size_t j = 0; j < n_edge; ++j) {
    const SegmentGas& sg = edgeGas(j);
    BoundaryLayerEdge& e = edge[j];
    e.x = sg.xm;
    e.radius = sg.r;
    e.u = std::max(sg.gas.u, 1.0e-3);
    e.rho = sg.gas.gas.rho;
    e.T = sg.gas.gas.T;
    e.p = sg.gas.gas.p;
    e.mach = sg.gas.mach;
    e.t_adiabatic_wall = sg.t_aw;
    e.viscosity = transport.mixture(sg.gas.gas.moleFractions(), sg.gas.gas.T,
                                    sg.gas.gas.cp_frozen).viscosity;
  }

  // Wall temperature per jacket segment: prescribed, or iterated.
  std::vector<double> t_wall(static_cast<std::size_t>(N),
                             coupled ? 0.8 * spec.wall.max_temperature : prescribed_wall);
  auto setWall = [&]() {
    const double upstream = spec.upstream_wall_temperature > 0.0
                                ? spec.upstream_wall_temperature
                                : t_wall.front();
    for (std::size_t j = 0; j < n_edge; ++j) {
      double tw = j < n_pad ? upstream : t_wall[j - n_pad];
      // The energy integral needs a driving difference; a wall at or above the
      // adiabatic temperature has none, so the guess is capped just below it.
      tw = std::min(tw, edge[j].t_adiabatic_wall - 1.0);
      edge[j].t_wall = tw;
    }
  };

  // --- passes ------------------------------------------------------------
  // One pass is a complete coolant march.  More are needed when something
  // the march depends on is itself an output of it: the boundary layer's wall
  // temperature, and a film drawn from the jacket's outlet.  Bartz with no
  // such film is a single pass.
  const int max_passes = ((use_bl && coupled) || film_auto_temperature) ? 8 : 1;
  constexpr double kTolerance = 0.5;   // K, on the wall and on the film temperature
  double t_film = film ? (film_auto_temperature ? spec.inlet_temperature : spec.film_temperature)
                       : 0.0;
  double t_film_used = t_film;
  CoolingResult res;
  BoundaryLayerSolution bl;
  int pass = 0;
  bool settled = false;
  for (; pass < max_passes; ++pass) {
    if (film) {
      updateFilm(t_film);
      t_film_used = t_film;
    }
    if (use_bl) {
      setWall();
      bl = marchBoundaryLayer(edge, referenceAt);
      res = runPass([&](int i, double t_wg) {
        const std::size_t j = n_pad + static_cast<std::size_t>(i);
        const BoundaryLayerEdge& e = edge[j];
        const double ts = eckertReferenceTemperature(e.T, t_wg, e.t_adiabatic_wall);
        return spec.hot_gas_multiplier *
               boundaryLayerFilmCoefficient(e, bl.points[j].enthalpy_thickness, t_wg,
                                            referenceAt(j, ts));
      });
    } else {
      res = runPass([&](int i, double t_wg) {
        const SegmentGas& sg = seg[static_cast<std::size_t>(i)];
        return bartzFilmCoefficient(bref, sg.eps, sg.gas.mach, sg.gas.gas.gamma_s, t_wg);
      });
    }
    double change = 0.0;
    if (use_bl && coupled) {
      for (int i = 0; i < N; ++i) {
        const double tw = res.stations[static_cast<std::size_t>(i)].t_wall_hot;
        change = std::max(change, std::abs(tw - t_wall[static_cast<std::size_t>(i)]));
        t_wall[static_cast<std::size_t>(i)] = tw;
      }
    }
    if (film_auto_temperature) {
      change = std::max(change, std::abs(res.coolant_outlet_temperature - t_film));
      t_film = res.coolant_outlet_temperature;
    }
    if (change < kTolerance) { settled = true; break; }
  }
  res.wall_iterations = std::min(pass + 1, max_passes);
  if (max_passes > 1 && !settled) {
    std::ostringstream os;
    os << "cooling: the coupled passes (boundary layer and wall"
       << (film_auto_temperature ? ", film temperature" : "") << ") did not settle to "
       << kTolerance << " K in " << max_passes << " passes";
    res.warnings.push_back(os.str());
  }
  if (use_bl) {
    for (int i = 0; i < N; ++i) {
      const auto& pt = bl.points[n_pad + static_cast<std::size_t>(i)];
      auto& st = res.stations[static_cast<std::size_t>(i)];
      st.bl_momentum_thickness = pt.momentum_thickness;
      st.bl_enthalpy_thickness = pt.enthalpy_thickness;
      st.bl_re_theta = pt.re_theta;
    }
    for (const auto& w : bl.warnings) res.warnings.push_back(w);
  }
  if (film) {
    res.has_film = true;
    res.film_mass_flow = spec.film_mass_flow;
    res.film_temperature = t_film_used;
    res.film_coolant_velocity = film_info.coolant_velocity;
    res.film_velocity_ratio = film_info.velocity_ratio;
    res.film_length_half = film_info.length_half;
    res.film_length_fifth = film_info.length_fifth;
    std::ostringstream os;
    if (film_info.clamped) {
      os << "film cooling: the gas-to-film velocity ratio " << film_info.velocity_ratio
         << " is outside the " << kHatchPapellMinVelocityRatio << "-" << kHatchPapellMaxVelocityRatio
         << " Hatch and Papell fitted; their velocity function is held at the nearer limit, and "
            "against their own data the correlation is unreliable out there";
      res.warnings.push_back(os.str());
    }
    if (!film_info.reaches_fifth) {
      std::ostringstream o2;
      o2 << "film cooling: the film's effectiveness stays above 0.2 to the end of the jacket; "
            "its reported reach is the jacket length";
      res.warnings.push_back(o2.str());
    }
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
