// SPDX-License-Identifier: MIT
#pragma once
/// \file BoundaryLayer.hpp
/// \brief Integral turbulent boundary layer along a nozzle wall.
///
/// WHY
/// ---
/// Bartz's closed-form correlation fixes the boundary layer's thickness by
/// assumption: it scales everything on the throat diameter and so has no idea
/// how far the layer has already grown.  Measured nozzles disagree with it in
/// exactly that way.  Against JPL TR 32-415 (air, 18 in of cooled approach) it
/// over-predicts by a median 1.4-1.6x; against NASA TN D-2832 (a LOX/GH2
/// rocket with 14.5 in of chamber) it is right in the chamber and 1.7x high at
/// the throat.  A boundary layer that is marched from where it starts carries
/// that history, and the same march gives the displacement and momentum
/// thicknesses the performance losses need.
///
/// MODEL
/// -----
/// The method is the classical integral one (Bartz, "Turbulent boundary-layer
/// heat transfer from rapidly accelerating flow of rocket combustion gases and
/// of heated air", Advances in Heat Transfer 2, 1-108, 1965; Elliott, Bartz &
/// Silver, JPL TR 32-387, 1963).  Two ordinary differential equations are
/// marched along the wall from the injector face:
///
///   momentum   dtheta/dx = Cf/2 - theta [ (2 + H - Me^2) dln(ue)/dx + dln(r)/dx ]
///   energy     dDelta/dx = St   - Delta [ dln(Ge dH)/dx + dln(r)/dx ]
///
/// with theta the momentum thickness, Delta the enthalpy (energy) thickness,
/// Ge = rho_e u_e the edge mass flux and dH = cp* (T_aw - T_w) the driving
/// enthalpy difference.  The edge density gradient in the momentum equation is
/// eliminated with the isentropic relation dln(rho_e) = -Me^2 dln(u_e).
///
/// CLOSURE (fixed before any comparison with data, and not tuned to it)
/// -------------------------------------------------------------------
///   * velocity profile u/u_e = (y/delta)^(1/7); temperature by Crocco-Busemann,
///     T = T_w + (T_aw - T_w)(u/u_e) - (T_aw - T_e)(u/u_e)^2.  The shape factor
///     H = delta*/theta follows by quadrature of that pair.
///   * properties at Eckert's reference temperature
///     T* = T_e + 0.5 (T_w - T_e) + 0.22 (T_aw - T_e), at the local static
///     pressure and edge composition (the same reference state NASA TN D-2832
///     correlated its data at).
///   * skin friction, the 1/7-law (Blasius) flat-plate relation with reference
///     properties:  Cf/2 = 0.0128 (rho*/rho_e) (G_e theta / mu*)^-1/4
///   * Stanton number by the Colburn analogy applied to the energy thickness:
///     St = 0.0128 (rho*/rho_e) (G_e Delta / mu*)^-1/4 Pr*^-2/3
///   * film coefficient h_g = St G_e cp*, per unit (T_aw - T_w).
///
/// On a flat plate the march reproduces the textbook 1/7-law result
/// theta/x = 0.036 Re_x^-0.2 exactly; that is checked in the unit tests.
///
/// WHAT IT DOES NOT MODEL
/// ----------------------
/// The closure is turbulent everywhere.  Laminarisation by strong acceleration
/// (acceleration parameter K = (nu/u^2) du/dx above about 3e-6) and laminar or
/// transitional layers (momentum-thickness Reynolds number below a few
/// hundred) are reported as warnings, not modelled.  JPL's two low-pressure
/// tests sit in that regime and this model over-predicts them as badly as
/// Bartz does.  No injector-face recirculation, film cooling (that is
/// HeatTransfer's film model), roughness, or wall curvature terms.

#include <functional>
#include <string>
#include <vector>

namespace ignis {

/// Free-stream (core-flow) and wall conditions at one axial station.
struct BoundaryLayerEdge {
  double x = 0.0;                 ///< m, along the axis from the injector face
  double radius = 0.0;            ///< m, wall radius
  double u = 0.0;                 ///< m/s, edge velocity
  double rho = 0.0;               ///< kg/m^3, edge density
  double T = 0.0;                 ///< K, edge static temperature
  double p = 0.0;                 ///< Pa, edge static pressure
  double mach = 0.0;
  double viscosity = 0.0;         ///< Pa s at the edge state (for the acceleration parameter)
  double t_adiabatic_wall = 0.0;  ///< K
  double t_wall = 0.0;            ///< K
};

/// Gas properties at a station's reference temperature.
struct ReferenceState {
  double rho = 0.0;       ///< kg/m^3
  double viscosity = 0.0; ///< Pa s
  double cp = 0.0;        ///< J/(kg K)
  double prandtl = 0.0;
};

/// Properties of the gas at station `index` (its pressure and composition) and
/// temperature `T`.  Supplied by the caller so that this module needs no
/// knowledge of which mixture or transport model is in use.
using ReferenceProperties = std::function<ReferenceState(std::size_t index, double T)>;

/// Eckert's reference temperature, K.
double eckertReferenceTemperature(double t_edge, double t_wall, double t_adiabatic_wall);

/// Shape factor H = delta*/theta of a 1/n power-law velocity profile with a
/// Crocco-Busemann temperature profile.
/// \param tw_te   T_wall / T_edge
/// \param taw_te  T_adiabatic_wall / T_edge
double compressibleShapeFactor(double tw_te, double taw_te, double exponent = 7.0);

/// One converged station of the march.
struct BoundaryLayerPoint {
  double x = 0.0;
  double momentum_thickness = 0.0;      ///< theta, m
  double enthalpy_thickness = 0.0;      ///< Delta, m
  double displacement_thickness = 0.0;  ///< delta* = H theta, m
  double shape_factor = 0.0;            ///< H
  double skin_friction = 0.0;           ///< Cf (twice the half value), on rho_e u_e^2 / 2
  double stanton = 0.0;                 ///< on G_e cp*
  double h_gas = 0.0;                   ///< W/(m^2 K), per (T_aw - T_w)
  double t_reference = 0.0;             ///< K
  double re_theta = 0.0;                ///< G_e theta / mu*
  double acceleration = 0.0;            ///< K = (mu_e / (rho_e u_e^2)) du_e/dx
};

/// Options for the march.
///
/// STARTING THICKNESS.  The layer is turbulent from its first station, as in
/// Elliott, Bartz & Silver (growth "from zero thickness") and Wang & Luong
/// (J. Thermophysics 8(3), 1994: "the turbulent boundary layer started from the
/// injector faceplate").  It cannot literally start at zero -- the closure's
/// Cf and St go as thickness^-1/4, so a zero start puts an infinite heat flux
/// at the first station, and a "small" start puts whatever the arbitrary value
/// implies there.  By default the march therefore starts from the thinnest
/// layer that can be turbulent at all: a momentum-thickness Reynolds number of
/// 320, Preston's minimum (J. Fluid Mech. 3, 373-384, 1958), with the energy
/// thickness equal to it (the "smooth start" of JPL TR 32-387, zeta_0 = 1).
/// That is a physical bound fixed before any comparison, not a fitted value;
/// for the engines in this repository it is a few hundredths of a millimetre
/// and the layer forgets it within millimetres.
struct BoundaryLayerOptions {
  /// Starting thicknesses, m.  Zero (the default) means "the minimum turbulent
  /// layer" described above; a positive value is used as given.
  double initial_momentum_thickness = 0.0;
  double initial_enthalpy_thickness = 0.0;
  double start_re_theta = 320.0;  ///< Re_theta of the default starting layer
  int substeps = 4;               ///< RK4 steps per station interval
  double profile_exponent = 7.0;  ///< n in u/u_e = (y/delta)^(1/n)
  /// Validity flags: above / below these the turbulent closure is suspect.
  double laminarisation_acceleration = 3.0e-6;
  double minimum_turbulent_re_theta = 320.0;   ///< Preston (1958)
};

/// Result of a march.
struct BoundaryLayerSolution {
  std::vector<BoundaryLayerPoint> points;   ///< one per edge station
  double max_acceleration = 0.0;
  double max_acceleration_x = 0.0;
  double min_re_theta = 0.0;                ///< smallest Re_theta of any station
  std::vector<std::string> warnings;

  /// Thicknesses interpolated linearly in x (clamped at the ends).
  BoundaryLayerPoint at(double x) const;
};

/// March the momentum and energy integrals along `edge` (strictly increasing x).
BoundaryLayerSolution marchBoundaryLayer(const std::vector<BoundaryLayerEdge>& edge,
                                         const ReferenceProperties& properties,
                                         const BoundaryLayerOptions& options = {});

/// Film coefficient of the closure at one station, for a given enthalpy
/// thickness and wall temperature, W/(m^2 K).  This is what a coupled wall
/// balance calls: the thickness comes from the march, the wall temperature is
/// the trial value.
double boundaryLayerFilmCoefficient(const BoundaryLayerEdge& edge, double enthalpy_thickness,
                                    double t_wall, const ReferenceState& reference);

}  // namespace ignis
