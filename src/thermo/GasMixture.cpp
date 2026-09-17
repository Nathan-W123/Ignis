// SPDX-License-Identifier: MIT
#include "ignis/thermo/GasMixture.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>

namespace ignis {

std::string toString(CompositionModel m) {
  return m == CompositionModel::kFrozen ? "frozen" : "equilibrium";
}

CompositionModel compositionModelFromString(const std::string& s) {
  if (s == "frozen") return CompositionModel::kFrozen;
  if (s == "equilibrium" || s == "shifting") return CompositionModel::kEquilibrium;
  throw ConfigError("composition model must be 'frozen' or 'equilibrium', got '" + s + "'");
}

double GasMixture::massResidual(const Eigen::VectorXd& n) const {
  return n.dot(db_->molarMasses()) - 1.0;
}

Eigen::VectorXd GasMixture::normalize(const Eigen::VectorXd& n) const {
  const double m = n.dot(db_->molarMasses());
  if (!(m > 0.0)) throw RangeError("cannot normalise a mixture with non-positive mass");
  return n / m;
}

Eigen::VectorXd GasMixture::fromMoleFractions(const Eigen::VectorXd& X) const {
  if (X.size() != static_cast<Eigen::Index>(db_->size()))
    throw ConfigError("mole-fraction vector has the wrong length");
  if ((X.array() < -1e-12).any()) throw RangeError("negative mole fraction supplied");
  const double sum = X.sum();
  if (!(sum > 0.0)) throw RangeError("mole fractions sum to zero");
  const double M = (X.dot(db_->molarMasses())) / sum;  // kg/mol
  return (X / sum) / M;                                 // mol per kg
}

Eigen::VectorXd GasMixture::fromMassFractions(const Eigen::VectorXd& Y) const {
  if (Y.size() != static_cast<Eigen::Index>(db_->size()))
    throw ConfigError("mass-fraction vector has the wrong length");
  if ((Y.array() < -1e-12).any()) throw RangeError("negative mass fraction supplied");
  const double sum = Y.sum();
  if (!(sum > 0.0)) throw RangeError("mass fractions sum to zero");
  return (Y / sum).cwiseQuotient(db_->molarMasses());
}

double GasMixture::molarMass(const Eigen::VectorXd& n) const {
  const double nt = n.sum();
  if (!(nt > 0.0)) throw RangeError("mixture has no moles");
  return n.dot(db_->molarMasses()) / nt;
}

double GasMixture::enthalpy(const Eigen::VectorXd& n, double T) const {
  double h = 0.0;
  for (std::size_t j = 0; j < db_->size(); ++j) {
    const double nj = n(static_cast<Eigen::Index>(j));
    if (nj != 0.0) h += nj * (*db_)[j].h(T);
  }
  return h;
}

double GasMixture::internalEnergy(const Eigen::VectorXd& n, double T) const {
  return enthalpy(n, T) - n.sum() * constants::R_universal * T;
}

double GasMixture::cpFrozen(const Eigen::VectorXd& n, double T) const {
  double cp = 0.0;
  for (std::size_t j = 0; j < db_->size(); ++j) {
    const double nj = n(static_cast<Eigen::Index>(j));
    if (nj != 0.0) cp += nj * (*db_)[j].cp(T);
  }
  return cp;
}

double GasMixture::entropy(const Eigen::VectorXd& n, double T, double p) const {
  if (!(p > 0.0)) throw RangeError("entropy requires positive pressure");
  const double nt = n.sum();
  if (!(nt > 0.0)) throw RangeError("entropy requires a non-empty mixture");
  const double lnp = std::log(p / constants::p_reference);
  double s = 0.0;
  for (std::size_t j = 0; j < db_->size(); ++j) {
    const double nj = n(static_cast<Eigen::Index>(j));
    if (nj <= 0.0) continue;  // the limit n ln n -> 0 makes absent species contribute nothing
    const double Xj = nj / nt;
    s += nj * ((*db_)[j].s0(T) - constants::R_universal * (std::log(Xj) + lnp));
  }
  return s;
}

Eigen::VectorXd GasMixture::elementMoles(const Eigen::VectorXd& n) const {
  return db_->elementMatrix() * n;
}

GasState GasMixture::frozenState(const Eigen::VectorXd& n, double T, double p) const {
  requireInRange(T);
  if (!(p > 0.0)) throw RangeError("frozenState requires positive pressure");
  GasState st;
  st.model = CompositionModel::kFrozen;
  st.T = T;
  st.p = p;
  st.n = n;
  st.n_total = n.sum();
  if (!(st.n_total > 0.0)) throw RangeError("frozenState requires a non-empty mixture");
  st.M = molarMass(n);
  st.R = constants::R_universal / st.M;
  st.rho = p / (st.R * T);
  st.v = 1.0 / st.rho;
  st.h = enthalpy(n, T);
  st.u = st.h - p * st.v;
  st.s = entropy(n, T, p);
  st.g = st.h - T * st.s;
  st.cp_frozen = cpFrozen(n, T);
  st.cv_frozen = st.cp_frozen - st.n_total * constants::R_universal;
  if (!(st.cv_frozen > 0.0))
    throw RangeError("non-positive frozen cv -- species data inconsistent");
  st.gamma_frozen = st.cp_frozen / st.cv_frozen;
  st.cp_eff = st.cp_frozen;
  st.cv_eff = st.cv_frozen;
  st.gamma_eff = st.gamma_frozen;
  st.dlnV_dlnT_p = 1.0;
  st.dlnV_dlnp_T = -1.0;
  st.gamma_s = st.gamma_frozen;
  st.a = std::sqrt(st.gamma_s * p / st.rho);
  return st;
}

void GasMixture::requireInRange(double T) const {
  for (const auto& s : db_->species()) {
    if (!s.inRange(T)) {
      std::ostringstream os;
      os << "temperature " << T << " K outside the fit range of species " << s.name() << " ["
         << s.tMin() << ", " << s.tMax() << "] K";
      throw RangeError(os.str());
    }
  }
}

namespace {
/// Safeguarded Newton with bisection fallback for monotone f(T) = 0.
double solveMonotone(const std::function<double(double)>& f,
                     const std::function<double(double)>& dfdT, double T_guess, double T_lo,
                     double T_hi, double tol, const char* what) {
  double lo = T_lo, hi = T_hi;
  double f_lo = f(lo), f_hi = f(hi);
  if (f_lo * f_hi > 0.0) {
    std::ostringstream os;
    os << what << ": target not bracketed on [" << lo << ", " << hi << "] K (residuals "
       << f_lo << ", " << f_hi << ")";
    throw ConvergenceError(os.str());
  }
  double T = std::min(std::max(T_guess, lo), hi);
  for (int it = 0; it < 200; ++it) {
    const double fT = f(T);
    if (std::abs(fT) < tol) return T;
    if (fT * f_lo > 0.0) { lo = T; f_lo = fT; } else { hi = T; f_hi = fT; }
    const double d = dfdT(T);
    double T_new = (d != 0.0) ? T - fT / d : 0.5 * (lo + hi);
    if (!(T_new > lo && T_new < hi) || !std::isfinite(T_new)) T_new = 0.5 * (lo + hi);
    if (std::abs(T_new - T) < 1e-12 * std::abs(T)) return T_new;
    T = T_new;
  }
  std::ostringstream os;
  os << what << ": failed to converge within 200 iterations (bracket width "
     << (hi - lo) << " K)";
  throw ConvergenceError(os.str());
}
}  // namespace

double GasMixture::temperatureFromEnthalpy(const Eigen::VectorXd& n, double h,
                                           double T_guess) const {
  const double lo = db_->tMinCommon(), hi = db_->tMaxCommon();
  const double scale = std::max(1.0, std::abs(h));
  return solveMonotone([&](double T) { return enthalpy(n, T) - h; },
                       [&](double T) { return cpFrozen(n, T); }, T_guess, lo, hi,
                       1e-10 * scale, "temperatureFromEnthalpy");
}

double GasMixture::temperatureFromInternalEnergy(const Eigen::VectorXd& n, double u,
                                                 double T_guess) const {
  const double lo = db_->tMinCommon(), hi = db_->tMaxCommon();
  const double nR = n.sum() * constants::R_universal;
  const double scale = std::max(1.0, std::abs(u));
  return solveMonotone([&](double T) { return internalEnergy(n, T) - u; },
                       [&](double T) { return cpFrozen(n, T) - nR; }, T_guess, lo, hi,
                       1e-10 * scale, "temperatureFromInternalEnergy");
}

double GasMixture::temperatureFromEntropy(const Eigen::VectorXd& n, double s, double p,
                                          double T_guess) const {
  const double lo = db_->tMinCommon(), hi = db_->tMaxCommon();
  const double scale = std::max(1.0, std::abs(s));
  return solveMonotone([&](double T) { return entropy(n, T, p) - s; },
                       [&](double T) { return cpFrozen(n, T) / T; }, T_guess, lo, hi,
                       1e-11 * scale, "temperatureFromEntropy");
}

}  // namespace ignis
