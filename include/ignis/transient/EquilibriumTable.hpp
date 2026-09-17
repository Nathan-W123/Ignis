// SPDX-License-Identifier: MIT
#pragma once
/// \file EquilibriumTable.hpp
/// \brief Pre-tabulated equilibrium chamber properties for transient work.
///
/// A zero-dimensional transient needs the chamber equation of state thousands
/// of times per run, and a full Gibbs minimisation at every Runge-Kutta stage
/// would dominate the cost.  Ignis therefore tabulates the equilibrium
/// properties of the product mixture on a (mixture ratio, temperature,
/// ln pressure) grid and interpolates trilinearly.
///
/// The table stores, for each node:
///   * specific internal energy u  (J/kg, absolute)
///   * mean molar mass M           (kg/mol)
///   * isentropic exponent gamma_s
///   * equilibrium cv              (J/(kg K))
/// plus, per mixture ratio, the reactant enthalpy and the heat of combustion
/// referred to 298.15 K.
///
/// INTERPOLATION ERROR IS MEASURED, NOT ASSUMED.  `measureError()` samples
/// random points inside the grid, solves them exactly, and reports the error
/// distribution; the transient application prints it and the documentation
/// quotes it.

#include <string>
#include <vector>

#include "ignis/combustion/Propellant.hpp"
#include "ignis/equilibrium/Equilibrium.hpp"

namespace ignis {

/// Table axis definition.
struct TableGrid {
  double mr_min = 1.0, mr_max = 6.0;
  int mr_points = 17;
  double t_min = 300.0, t_max = 4200.0;
  int t_points = 49;
  double p_min = 1.0e4, p_max = 2.0e7;
  int p_points = 13;
};

/// Measured accuracy of the interpolation.
struct InterpolationError {
  int samples = 0;
  double max_rel_temperature = 0.0;
  double rms_rel_temperature = 0.0;
  double max_rel_pressure = 0.0;
  double rms_rel_pressure = 0.0;
  double max_rel_molar_mass = 0.0;
  double max_rel_cstar = 0.0;
  std::string summary() const;
};

/// Trilinear equilibrium property table.
class EquilibriumTable {
 public:
  /// Build the table.  `base` supplies the propellant pair and inlet
  /// temperatures; only its mixture ratio is swept.
  static EquilibriumTable build(const EquilibriumSolver& solver, const PropellantMixture& base,
                                const TableGrid& grid, bool verbose = false);

  /// Write the table to CSV for inspection and plotting.  This is an export
  /// format, not a reload format: rebuilding takes a few seconds and keeps the
  /// table guaranteed consistent with the species data in use.
  void exportCsv(const std::string& path) const;

  const TableGrid& grid() const { return grid_; }

  /// Specific internal energy of the equilibrium mixture, J/kg.
  double internalEnergy(double mr, double T, double p) const;
  /// Mean molar mass, kg/mol.
  double molarMass(double mr, double T, double p) const;
  /// Isentropic exponent.
  double gammaS(double mr, double T, double p) const;
  /// Equilibrium specific heat at constant volume, J/(kg K).
  double cv(double mr, double T, double p) const;
  /// Ideal characteristic velocity at the adiabatic flame condition, m/s.
  /// This is the full variable-property sonic solution, not a constant-gamma
  /// estimate.
  double cStar(double mr, double p) const;

  /// Characteristic velocity of a chamber that is *not* at its adiabatic flame
  /// temperature, m/s.
  ///
  /// The constant-gamma choked-flow result
  ///     c*_cg = sqrt(R T / gamma) ((gamma+1)/2)^((gamma+1)/(2(gamma-1)))
  /// is corrected by the factor kappa(O/F, p) = c*_exact / c*_cg measured at
  /// the adiabatic flame condition, so the value is exact there and degrades
  /// gracefully as the chamber cools during a start or a shutdown.
  double cStarAt(double mr, double T, double p) const;

  /// The ratio c*_exact / c*_constant-gamma at the adiabatic flame condition.
  double cStarCorrection(double mr, double p) const;
  /// Adiabatic flame temperature, K.
  double flameTemperature(double mr, double p) const;
  /// Reactant specific enthalpy at the stored inlet temperatures, J/kg.
  double reactantEnthalpy(double mr) const;
  /// Heat of combustion referred to 298.15 K and 1 bar, J/kg (positive).
  double heatOfCombustion(double mr) const;

  /// Solve the chamber equation of state: given the mixture ratio, specific
  /// internal energy and density, return the temperature and pressure.
  /// \throws ConvergenceError if the state lies outside the table.
  void solveState(double mr, double u, double rho, double& T, double& p) const;

  /// Sample random interior points and compare against exact solves.
  InterpolationError measureError(const EquilibriumSolver& solver, const PropellantMixture& base,
                                  int samples, unsigned seed) const;

  double mrMin() const { return grid_.mr_min; }
  double mrMax() const { return grid_.mr_max; }
  double tMin() const { return grid_.t_min; }
  double tMax() const { return grid_.t_max; }

 private:
  std::size_t idx(int i, int j, int k) const {
    return (static_cast<std::size_t>(i) * static_cast<std::size_t>(grid_.t_points) +
            static_cast<std::size_t>(j)) * static_cast<std::size_t>(grid_.p_points) +
           static_cast<std::size_t>(k);
  }
  double interp3(const std::vector<double>& f, double mr, double T, double p) const;
  double interp2(const std::vector<double>& f, double mr, double p) const;
  double interp1mr(const std::vector<double>& f, double mr) const;

  TableGrid grid_;
  std::vector<double> mr_, t_, p_, lnp_;
  std::vector<double> u_, molar_, gamma_, cv_;         // (mr, T, p)
  std::vector<double> cstar_, tflame_, kappa_;         // (mr, p)
  std::vector<double> h_react_, q_comb_;               // (mr)
};

}  // namespace ignis
