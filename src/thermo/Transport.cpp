// SPDX-License-Identifier: MIT
#include "ignis/thermo/Transport.hpp"

#include <cmath>
#include <sstream>

namespace ignis {
namespace {
constexpr double kDebye = 3.33564e-30;      // C m per Debye
constexpr double kFourPiEps0 = 1.11265005545e-10;  // C^2 / (J m)
constexpr double kAngstrom = 1.0e-10;       // m
}  // namespace

double omega22LJ(double t_star) {
  // Neufeld, Janzen & Aziz (1972), 3-term form.
  return 1.16145 * std::pow(t_star, -0.14874) + 0.52487 * std::exp(-0.77320 * t_star) +
         2.16178 * std::exp(-2.43787 * t_star);
}

double omega11LJ(double t_star) {
  return 1.06036 * std::pow(t_star, -0.15610) + 0.19300 * std::exp(-0.47635 * t_star) +
         1.03587 * std::exp(-1.52996 * t_star) + 1.76474 * std::exp(-3.89411 * t_star);
}

TransportModel::TransportModel(const SpeciesDatabase& db) : db_(&db) {
  delta_star_.assign(db.size(), 0.0);
  has_data_.assign(db.size(), false);
  t_min_valid_ = 0.0;
  t_max_valid_ = 1.0e30;
  for (std::size_t j = 0; j < db.size(); ++j) {
    const auto& tr = db[j].transport();
    has_data_[j] = tr.valid;
    if (!tr.valid) continue;
    const double eps = tr.well_depth * constants::k_B;          // J
    const double sigma = tr.diameter * kAngstrom;               // m
    const double mu_d = tr.dipole * kDebye;                     // C m
    delta_star_[j] = 0.5 * mu_d * mu_d / (kFourPiEps0 * eps * sigma * sigma * sigma);
    // Neufeld correlation is fitted for 0.3 <= T* <= 100.
    t_min_valid_ = std::max(t_min_valid_, 0.3 * tr.well_depth);
    t_max_valid_ = std::min(t_max_valid_, 100.0 * tr.well_depth);
  }
}

double TransportModel::speciesViscosity(std::size_t j, double T) const {
  if (!has_data_[j])
    throw RangeError("no Lennard-Jones data for species " + (*db_)[j].name());
  const auto& tr = (*db_)[j].transport();
  const double t_star = T / tr.well_depth;
  const double omega = omega22LJ(t_star) + 0.2 * delta_star_[j] * delta_star_[j] / t_star;
  const double m = (*db_)[j].molarMass() / constants::N_A;   // kg per molecule
  const double sigma = tr.diameter * kAngstrom;
  return (5.0 / 16.0) * std::sqrt(constants::pi * m * constants::k_B * T) /
         (constants::pi * sigma * sigma * omega);
}

double TransportModel::speciesConductivity(std::size_t j, double T) const {
  const double mu = speciesViscosity(j, T);
  const double cp_molar = (*db_)[j].cp(T);
  const double cv_molar = cp_molar - constants::R_universal;
  // Modified Eucken; exact monatomic limit when cv = 3R/2.
  return (mu / (*db_)[j].molarMass()) * (1.32 * cv_molar + 1.77 * constants::R_universal);
}

TransportResult TransportModel::mixture(const Eigen::VectorXd& X, double T,
                                        double cp_mass) const {
  const Eigen::Index N = static_cast<Eigen::Index>(db_->size());
  if (X.size() != N) throw ConfigError("transport: mole-fraction vector has the wrong length");
  if (!(T > 0.0)) throw RangeError("transport: temperature must be positive");

  std::vector<Eigen::Index> idx;
  double covered = 0.0, total = X.sum();
  if (!(total > 0.0)) throw RangeError("transport: mole fractions sum to zero");
  for (Eigen::Index j = 0; j < N; ++j) {
    if (X(j) > 0.0 && has_data_[static_cast<std::size_t>(j)]) {
      idx.push_back(j);
      covered += X(j);
    }
  }
  if (idx.empty()) throw RangeError("transport: no species in the mixture have LJ data");

  const std::size_t n = idx.size();
  std::vector<double> x(n), mu(n), lam(n), M(n);
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t j = static_cast<std::size_t>(idx[k]);
    x[k] = X(idx[k]) / covered;                 // renormalised over covered species
    mu[k] = speciesViscosity(j, T);
    lam[k] = speciesConductivity(j, T);
    M[k] = (*db_)[j].molarMass();
  }

  TransportResult out;
  double mu_mix = 0.0, lam_mix = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    double denom = 0.0;
    for (std::size_t l = 0; l < n; ++l) {
      const double r = std::sqrt(mu[k] / mu[l]) * std::pow(M[l] / M[k], 0.25);
      const double phi = (1.0 + r) * (1.0 + r) / std::sqrt(8.0 * (1.0 + M[k] / M[l]));
      denom += x[l] * phi;
    }
    mu_mix += x[k] * mu[k] / denom;
    lam_mix += x[k] * lam[k] / denom;
  }
  out.viscosity = mu_mix;
  out.conductivity = lam_mix;
  out.covered_mole_fraction = covered / total;
  out.prandtl = (lam_mix > 0.0) ? mu_mix * cp_mass / lam_mix : 0.0;
  return out;
}

}  // namespace ignis
