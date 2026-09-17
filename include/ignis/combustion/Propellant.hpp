// SPDX-License-Identifier: MIT
#pragma once
/// \file Propellant.hpp
/// \brief Reactant definitions and propellant mixing.
///
/// A propellant carries the *absolute* molar enthalpy of its storage state, on
/// the same scale as the NASA-7 species polynomials.  The reactant enthalpy of
/// a mixture is therefore simply the mole-weighted sum, and feeding it to the
/// constant-enthalpy equilibrium problem gives the adiabatic flame temperature
/// with no separate heat-of-reaction term.

#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis {

/// One storable reactant.
struct Propellant {
  std::string name;
  std::string description;
  std::string role;                 ///< "oxidizer", "fuel" or "diluent"
  bool liquid = true;
  std::map<std::string, int> composition;
  double molar_mass = 0.0;          ///< kg/mol
  double reference_temperature = 0.0;  ///< K
  /// Absolute molar enthalpy at reference_temperature, J/mol.  Only used for
  /// liquids; gases evaluate the NASA-7 polynomial of `species` instead.
  double reference_enthalpy = 0.0;
  double cp_liquid = 0.0;           ///< J/(mol K), for the sensible correction
  double density = 0.0;             ///< kg/m^3 at reference_temperature
  std::string species;              ///< gas-phase species name (gases only)
  std::string coolant_table;        ///< real-fluid table name, if any
  /// Largest |T_in - T_ref| for which the linear sensible correction is
  /// accepted, K.  Beyond this the request is rejected rather than silently
  /// extrapolated.
  double max_temperature_offset = 20.0;
  std::string source;

  /// Absolute molar enthalpy at an inlet temperature, J/mol.
  double molarEnthalpy(double T_inlet, const SpeciesDatabase& db) const;
};

/// The shipped propellant database.
class PropellantLibrary {
 public:
  static PropellantLibrary loadYaml(const std::string& path);
  const Propellant& at(const std::string& name) const;
  bool has(const std::string& name) const { return by_name_.count(name) > 0; }
  std::vector<std::string> names() const;
  const std::string& path() const { return path_; }

 private:
  std::map<std::string, Propellant> by_name_;
  std::string path_;
};

/// Locate the shipped propellant library (same search rules as the species
/// database).
std::string findDefaultPropellantLibrary();

/// A fuel/oxidiser pair at a commanded mixture ratio and inlet temperatures.
class PropellantMixture {
 public:
  PropellantMixture(Propellant oxidizer, Propellant fuel, double mixture_ratio,
                    double T_oxidizer, double T_fuel);

  const Propellant& oxidizer() const { return ox_; }
  const Propellant& fuel() const { return fuel_; }
  /// Oxidiser-to-fuel mass ratio.
  double mixtureRatio() const { return mr_; }
  double oxidizerTemperature() const { return T_ox_; }
  double fuelTemperature() const { return T_fuel_; }
  /// Mass fraction of oxidiser in the total flow.
  double oxidizerMassFraction() const { return mr_ / (1.0 + mr_); }
  double fuelMassFraction() const { return 1.0 / (1.0 + mr_); }

  /// Element mole numbers per kg of propellant mixture, in the row order of
  /// `db.elements()`.
  Eigen::VectorXd elementMoles(const SpeciesDatabase& db) const;
  /// Specific enthalpy of the unburned reactants, J/kg.
  double enthalpy(const SpeciesDatabase& db) const;
  /// Mixture ratio at which all carbon burns to CO2 and all hydrogen to H2O.
  double stoichiometricMixtureRatio() const;
  /// phi = (F/O) / (F/O)_stoichiometric.  phi > 1 is fuel rich.
  double equivalenceRatio() const { return stoichiometricMixtureRatio() / mr_; }
  /// Bulk density of the propellant combination, kg/m^3 (volume-weighted).
  double bulkDensity() const;

  /// Return a copy at a different mixture ratio.
  PropellantMixture withMixtureRatio(double mr) const;

 private:
  Propellant ox_, fuel_;
  double mr_;
  double T_ox_, T_fuel_;
};

}  // namespace ignis
