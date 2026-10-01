// SPDX-License-Identifier: MIT
#include "ignis/thermal/BoundaryLayerLosses.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/thermal/HeatTransfer.hpp"

namespace ignis {

BoundaryLayerLossResult computeBoundaryLayerLosses(
    const NozzleFlow& flow, const NozzleGeometry& geom, const ChamberResult& chamber,
    const TransportModel& transport, const std::function<double(double)>& wall_temperature,
    const BoundaryLayerLossOptions& opts) {
  if (!(opts.uncooled_wall_temperature > 0.0))
    throw ConfigError("boundary-layer losses: the uncooled wall temperature must be positive");
  if (opts.stations < 20) throw ConfigError("boundary-layer losses: at least 20 stations");

  // Recovery factor Pr^(1/3) at the chamber Prandtl number, as in the
  // cooling solve, so both define the driving temperature the same way.
  const auto tr0 = transport.mixture(chamber.state.moleFractions(), chamber.state.T,
                                     chamber.state.cp_frozen);
  const double pr0 = tr0.prandtl;
  const GasMixture& mix = flow.mixture();

  const int n = opts.stations;
  const double x_end = geom.exitPosition();
  std::vector<BoundaryLayerEdge> edge(static_cast<std::size_t>(n));
  std::vector<ExpansionState> state(static_cast<std::size_t>(n));
  int cooled = 0;
  for (int i = 0; i < n; ++i) {
    BoundaryLayerEdge& e = edge[static_cast<std::size_t>(i)];
    e.x = x_end * i / (n - 1);
    e.radius = geom.radius(e.x);
    const double eps = std::max(geom.area(e.x) / geom.throatArea(), 1.0);
    const ExpansionState s = flow.atAreaRatio(eps, e.x > geom.throatPosition());
    state[static_cast<std::size_t>(i)] = s;
    e.u = std::max(s.u, 1.0e-3);
    e.rho = s.gas.rho;
    e.T = s.gas.T;
    e.p = s.gas.p;
    e.mach = s.mach;
    e.viscosity = transport.mixture(s.gas.moleFractions(), s.gas.T, s.gas.cp_frozen).viscosity;
    e.t_adiabatic_wall = recoveryTemperature(s.gas.T, s.mach, s.gas.gamma_s, pr0);
    double tw = wall_temperature ? wall_temperature(e.x) : 0.0;
    if (tw > 0.0) {
      ++cooled;
    } else {
      tw = opts.uncooled_wall_temperature;
    }
    // The energy integral needs a positive driving difference.
    e.t_wall = std::min(tw, e.t_adiabatic_wall - 1.0);
  }
  auto properties = [&](std::size_t i, double T) {
    const ExpansionState& s = state[i];
    ReferenceState r;
    r.rho = s.gas.p * s.gas.M / (constants::R_universal * T);
    r.cp = mix.cpFrozen(s.gas.n, T);
    const auto t = transport.mixture(s.gas.moleFractions(), T, r.cp);
    r.viscosity = t.viscosity;
    r.prandtl = t.prandtl;
    return r;
  };

  BoundaryLayerLossResult out;
  out.layer = marchBoundaryLayer(edge, properties);
  out.cooled_fraction = static_cast<double>(cooled) / n;
  const auto throat = out.layer.at(geom.throatPosition());
  const auto& exit = out.layer.points.back();
  out.throat_displacement_thickness = throat.displacement_thickness;
  out.throat_momentum_thickness = throat.momentum_thickness;
  out.exit_displacement_thickness = exit.displacement_thickness;
  out.exit_momentum_thickness = exit.momentum_thickness;
  out.exit_re_theta = exit.re_theta;

  const double rt = geom.throatRadius(), re = geom.exitRadius();
  if (!(out.throat_displacement_thickness < rt && out.exit_displacement_thickness < re))
    throw ConvergenceError("boundary-layer losses: the displacement thickness reached the axis; "
                           "the thin-layer corrections do not apply");
  out.discharge_coefficient = std::pow(1.0 - out.throat_displacement_thickness / rt, 2);
  out.geometric_area_ratio = geom.expansionRatio();
  out.effective_area_ratio = std::pow((re - out.exit_displacement_thickness) /
                                      (rt - out.throat_displacement_thickness), 2);
  for (const auto& w : out.layer.warnings) out.warnings.push_back(w);
  if (cooled < n) {
    std::ostringstream os;
    os << "boundary-layer losses: " << n - cooled << " of " << n
       << " stations have no jacket wall temperature and assume "
       << opts.uncooled_wall_temperature << " K";
    out.warnings.push_back(os.str());
  }
  return out;
}

}  // namespace ignis
