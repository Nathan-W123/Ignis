// SPDX-License-Identifier: MIT
#include "ignis/combustion/Chamber.hpp"

#include <iomanip>
#include <sstream>

namespace ignis {

double ChamberResult::throatAreaFor(double mdot) const {
  if (!(mdot > 0.0)) throw ConfigError("chamber: mass flow must be positive");
  return mdot * c_star / state.p;
}

double ChamberResult::massFlowFor(double at) const {
  if (!(at > 0.0)) throw ConfigError("chamber: throat area must be positive");
  return state.p * at / c_star;
}

double ChamberResult::residenceTime(double chamber_volume, double mdot) const {
  if (!(chamber_volume > 0.0)) throw ConfigError("chamber: volume must be positive");
  if (!(mdot > 0.0)) throw ConfigError("chamber: mass flow must be positive");
  return chamber_volume * state.rho / mdot;
}

std::string ChamberResult::summary() const {
  std::ostringstream os;
  os << std::fixed;
  os << "chamber (" << toString(model) << " expansion basis)\n"
     << "  mixture ratio O/F     " << std::setprecision(4) << mixture_ratio
     << "   (stoichiometric " << stoichiometric_mixture_ratio << ", phi "
     << equivalence_ratio << ")\n"
     << "  reactant enthalpy     " << std::setprecision(1) << reactant_enthalpy * 1e-3
     << " kJ/kg\n"
     << "  chamber pressure      " << std::setprecision(4) << state.p * 1e-6 << " MPa\n"
     << "  flame temperature     " << std::setprecision(2) << state.T << " K\n"
     << "  molar mass            " << std::setprecision(4) << state.M * 1e3 << " g/mol\n"
     << "  gas constant          " << std::setprecision(2) << state.R << " J/(kg K)\n"
     << "  density               " << std::setprecision(4) << state.rho << " kg/m^3\n"
     << "  cp (frozen / shift)   " << std::setprecision(1) << state.cp_frozen << " / "
     << state.cp_eff << " J/(kg K)\n"
     << "  cv (frozen / shift)   " << state.cv_frozen << " / " << state.cv_eff << " J/(kg K)\n"
     << "  gamma (frozen)        " << std::setprecision(5) << state.gamma_frozen << "\n"
     << "  gamma_s (isentropic)  " << state.gamma_s << "\n"
     << "  speed of sound        " << std::setprecision(2) << state.a << " m/s\n"
     << "  entropy               " << std::setprecision(2) << state.s << " J/(kg K)\n"
     << "  c* ideal              " << std::setprecision(2) << c_star_ideal << " m/s\n"
     << "  eta_c*                " << std::setprecision(4) << eta_c_star << "\n"
     << "  c* actual             " << std::setprecision(2) << c_star << " m/s\n"
     << "  throat mass flux      " << std::setprecision(2) << throat_mass_flux
     << " kg/(m^2 s)\n"
     << "  throat temperature    " << throat_state.T << " K\n"
     << "  equilibrium residuals: element " << std::scientific << std::setprecision(2)
     << diagnostics.element_residual_rel << ", Gibbs " << diagnostics.gibbs_residual
     << ", enthalpy " << diagnostics.state_residual << ", mass " << diagnostics.mass_residual;
  return os.str();
}

CombustionChamber::CombustionChamber(const EquilibriumSolver& solver, CompositionModel model)
    : solver_(&solver), model_(model) {}

ChamberResult CombustionChamber::solve(const PropellantMixture& mix, double p_chamber,
                                       double eta_c_star) const {
  if (!(p_chamber > 0.0))
    throw ConfigError("chamber: pressure must be positive, got " + std::to_string(p_chamber));
  if (!(eta_c_star > 0.0 && eta_c_star <= 1.0))
    throw ConfigError("chamber: eta_c_star must lie in (0, 1], got " +
                      std::to_string(eta_c_star));

  ChamberResult res;
  res.model = model_;
  res.mixture_ratio = mix.mixtureRatio();
  res.stoichiometric_mixture_ratio = mix.stoichiometricMixtureRatio();
  res.equivalence_ratio = mix.equivalenceRatio();
  res.element_moles = mix.elementMoles(solver_->database());
  res.reactant_enthalpy = mix.enthalpy(solver_->database());

  const double element_mass = solver_->elementMass(res.element_moles);
  if (std::abs(element_mass - 1.0) > 1.0e-6) {
    std::ostringstream os;
    os << "chamber: the element vector represents " << element_mass
       << " kg instead of 1 kg; the propellant molar masses and the database atomic "
          "weights are inconsistent";
    throw ConfigError(os.str());
  }

  const auto eq = solver_->hp(res.element_moles, res.reactant_enthalpy, p_chamber);
  res.state = eq.state;
  res.diagnostics = eq.diagnostics;

  const NozzleFlow flow(*solver_, model_, ChamberReference{res.element_moles, res.state});
  res.c_star_ideal = flow.cStarIdeal();
  res.eta_c_star = eta_c_star;
  res.c_star = eta_c_star * res.c_star_ideal;
  res.throat_mass_flux = flow.throatMassFlux();
  res.throat_state = flow.throat().gas;
  res.throat_mach = flow.throat().mach;
  return res;
}

NozzleFlow CombustionChamber::makeFlow(const ChamberResult& res) const {
  return NozzleFlow(*solver_, model_, ChamberReference{res.element_moles, res.state});
}

}  // namespace ignis
