// SPDX-License-Identifier: MIT
#pragma once
/// \file BoundaryLayerLosses.hpp
/// \brief What the wall boundary layer costs the nozzle: discharge coefficient
///        and thrust deficit.
///
/// The quasi-1D core is inviscid.  The real flow carries a boundary layer whose
/// displacement thickness delta* narrows the passage the core sees and whose
/// momentum thickness theta is momentum the wall has taken out of the flow.
/// Marching the layer of BoundaryLayer.hpp from the injector face to the exit
/// plane gives both, and with them three corrections, the standard ones of the
/// boundary-layer method in JANNAF performance practice (and of Elliott, Bartz
/// & Silver, JPL TR 32-387, Section V):
///
///   throat   A*_eff = pi (r_t - delta*_t)^2,     C_d = A*_eff / A*
///            mdot   = C_d p_c A* / c*
///   exit     the core expands to A_e,eff / A*_eff, A_e,eff = pi (r_e - delta*_e)^2
///   thrust   F = lambda eta (mdot u_e - 2 pi r_e rho_e u_e^2 theta_e) + (p_e - p_a) A_e
///
/// with the edge state (u_e, rho_e, p_e) that of the core at the effective
/// area ratio.  The thin-layer form 2 pi r theta is used for the momentum
/// deficit; on the engines here the exit layer is 1-3 % of the exit radius.
///
/// A cold wall shrinks delta*: the gas next to it is denser than the edge gas
/// and carries more of the mass the velocity deficit removed.  Severe enough
/// cooling can make it NEGATIVE -- the layer then carries more mass than the
/// inviscid flow it replaced, C_d exceeds one and the effective exit is larger
/// than the geometric one (Elliott, Bartz & Silver's Fig. 11 shows this for
/// their sample nozzle).  Nothing here clips that; it is what the integral
/// definitions give.
///
/// WALL TEMPERATURE.  The layer needs one everywhere.  Inside a solved cooling
/// jacket it is the jacket's hot-wall temperature; elsewhere -- an uncooled
/// extension, or no cooling analysis at all -- it is `uncooled_wall_temperature`,
/// an assumption the caller should know they are making.

#include <functional>
#include <string>
#include <vector>

#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"
#include "ignis/thermal/BoundaryLayer.hpp"
#include "ignis/thermo/Transport.hpp"

namespace ignis {

struct BoundaryLayerLossOptions {
  /// Wall temperature where no jacket temperature is supplied, K.
  double uncooled_wall_temperature = 1000.0;
  int stations = 600;   ///< march stations, injector face to exit
};

struct BoundaryLayerLossResult {
  double discharge_coefficient = 1.0;
  double throat_displacement_thickness = 0.0;  ///< m
  double throat_momentum_thickness = 0.0;      ///< m
  double exit_displacement_thickness = 0.0;    ///< m
  double exit_momentum_thickness = 0.0;        ///< m
  double exit_re_theta = 0.0;
  double effective_area_ratio = 0.0;           ///< A_e,eff / A*_eff
  double geometric_area_ratio = 0.0;
  double cooled_fraction = 0.0;                ///< share of stations with a jacket wall temperature
  BoundaryLayerSolution layer;                 ///< the march itself, for plotting
  std::vector<std::string> warnings;
};

/// March the layer along the whole contour and reduce it to the corrections.
///
/// \param wall_temperature  T_w(x) where a jacket supplies it; return a
///                          non-positive value where it does not, and the
///                          uncooled temperature is used there.
BoundaryLayerLossResult computeBoundaryLayerLosses(
    const NozzleFlow& flow, const NozzleGeometry& geom, const ChamberResult& chamber,
    const TransportModel& transport, const std::function<double(double)>& wall_temperature,
    const BoundaryLayerLossOptions& options = {});

}  // namespace ignis
