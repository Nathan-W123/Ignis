// SPDX-License-Identifier: MIT
#include "ignis/thermal/CoolantFluid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/core/Version.hpp"

namespace ignis {

std::string toString(CoolantPhase p) {
  switch (p) {
    case CoolantPhase::kSupercriticalPressure: return "supercritical-pressure";
    case CoolantPhase::kLiquid: return "liquid";
    case CoolantPhase::kVapor: return "vapor";
    case CoolantPhase::kSupercriticalTemperature: return "supercritical-temperature";
  }
  return "?";
}

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  std::istringstream is(s);
  while (std::getline(is, cur, sep)) out.push_back(cur);
  return out;
}

double parseTagged(const std::string& line, const std::string& key) {
  const auto pos = line.find(key + "=");
  if (pos == std::string::npos)
    throw ConfigError("coolant table: missing '" + key + "' in header line");
  return std::atof(line.c_str() + pos + key.size() + 1);
}

}  // namespace

CoolantFluid CoolantFluid::loadCsv(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ConfigError("cannot open coolant table '" + path + "'");

  CoolantFluid f;
  {
    namespace fs = std::filesystem;
    f.name_ = fs::path(path).stem().string();
  }
  std::string line;
  std::vector<double> sat_p, sat_t;
  int nT = 0, np = 0;
  std::ostringstream prov;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] != '#') break;  // this is the column header line
    if (line.rfind("# source:", 0) == 0 || line.rfind("# generated", 0) == 0)
      prov << line.substr(2) << "; ";
    if (line.rfind("# critical:", 0) == 0) {
      f.t_crit_ = parseTagged(line, "Tc");
      f.p_crit_ = parseTagged(line, "pc");
      f.rho_crit_ = parseTagged(line, "rhoc");
      f.molar_mass_ = parseTagged(line, "M");
      f.acentric_ = parseTagged(line, "acentric");
    } else if (line.rfind("#sat ", 0) == 0) {
      for (const auto& tok : split(line.substr(5), ' ')) {
        const auto c = tok.find(':');
        if (c == std::string::npos) continue;
        sat_p.push_back(std::atof(tok.substr(0, c).c_str()));
        sat_t.push_back(std::atof(tok.substr(c + 1).c_str()));
      }
    } else if (line.rfind("#grid", 0) == 0) {
      nT = static_cast<int>(parseTagged(line, "nT"));
      np = static_cast<int>(parseTagged(line, "np"));
    }
  }
  if (nT <= 1 || np <= 1) throw ConfigError("coolant table '" + path + "': missing grid header");
  f.provenance_ = prov.str();

  const std::size_t n = static_cast<std::size_t>(nT) * static_cast<std::size_t>(np);
  f.rho_.resize(n); f.cp_.resize(n); f.h_.resize(n); f.mu_.resize(n); f.k_.resize(n);
  f.t_.resize(static_cast<std::size_t>(nT));
  f.p_.resize(static_cast<std::size_t>(np));

  std::size_t idx = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto tok = split(line, ',');
    if (tok.size() < 8) throw ConfigError("coolant table '" + path + "': malformed data row");
    if (idx >= n) throw ConfigError("coolant table '" + path + "': more rows than the grid");
    const int i = static_cast<int>(idx) / np;
    const int j = static_cast<int>(idx) % np;
    f.t_[static_cast<std::size_t>(i)] = std::atof(tok[0].c_str());
    f.p_[static_cast<std::size_t>(j)] = std::atof(tok[1].c_str());
    if (std::atoi(tok[7].c_str()) == 0) {
      // Preserve the hole; `at()` refuses to interpolate across it.
      f.rho_[idx] = f.cp_[idx] = f.h_[idx] = f.mu_[idx] = f.k_[idx] =
          std::numeric_limits<double>::quiet_NaN();
    } else {
      f.rho_[idx] = std::atof(tok[2].c_str());
      f.cp_[idx] = std::atof(tok[3].c_str());
      f.h_[idx] = std::atof(tok[4].c_str());
      f.mu_[idx] = std::atof(tok[5].c_str());
      f.k_[idx] = std::atof(tok[6].c_str());
    }
    ++idx;
  }
  if (idx != n)
    throw ConfigError("coolant table '" + path + "': expected " + std::to_string(n) +
                      " rows, found " + std::to_string(idx));

  f.lnp_.resize(f.p_.size());
  for (std::size_t j = 0; j < f.p_.size(); ++j) f.lnp_[j] = std::log(f.p_[j]);
  // Per-pressure-column valid temperature window.  A hole at the bottom of a
  // column is the melting line; one at the top would be the end of the fit.
  f.first_valid_.assign(f.p_.size(), -1);
  f.last_valid_.assign(f.p_.size(), -1);
  for (std::size_t j = 0; j < f.p_.size(); ++j) {
    for (std::size_t i = 0; i < f.t_.size(); ++i) {
      if (std::isfinite(f.rho_[i * f.p_.size() + j])) {
        if (f.first_valid_[j] < 0) f.first_valid_[j] = static_cast<int>(i);
        f.last_valid_[j] = static_cast<int>(i);
      }
    }
    if (f.first_valid_[j] < 0)
      throw ConfigError("coolant table '" + path + "': pressure column " + std::to_string(j) +
                        " has no valid data at all");
  }
  f.sat_t_.assign(f.p_.size(), -1.0);
  for (std::size_t j = 0; j < sat_p.size() && j < f.sat_t_.size(); ++j)
    f.sat_t_[j] = std::isfinite(sat_t[j]) ? sat_t[j] : -1.0;
  return f;
}

CoolantFluid CoolantFluid::load(const std::string& fluid_name) {
  namespace fs = std::filesystem;
  std::vector<std::string> tried;
  auto probe = [&](const fs::path& p) {
    tried.push_back(p.string());
    std::error_code ec;
    return fs::exists(p, ec);
  };
  const char* env = std::getenv("IGNIS_DATA_DIR");
  if (env != nullptr) {
    const fs::path p = fs::path(env) / "coolants" / (fluid_name + ".csv");
    if (probe(p)) return loadCsv(p.string());
  }
  {
    const fs::path p = fs::path(kDefaultDataDir) / "coolants" / (fluid_name + ".csv");
    if (probe(p)) return loadCsv(p.string());
  }
  for (const char* rel : {"data/coolants/", "../data/coolants/", "../../data/coolants/"}) {
    const fs::path p = fs::path(rel) / (fluid_name + ".csv");
    if (probe(p)) return loadCsv(p.string());
  }
  std::ostringstream os;
  os << "cannot locate the coolant property table for '" << fluid_name << "'. Tried:";
  for (const auto& t : tried) os << "\n  " << t;
  throw ConfigError(os.str());
}

double CoolantFluid::interp(const std::vector<double>& field, double T, double p) const {
  const std::size_t np = p_.size();
  auto it = std::upper_bound(t_.begin(), t_.end(), T);
  std::size_t i = static_cast<std::size_t>(std::distance(t_.begin(), it));
  i = (i == 0) ? 0 : i - 1;
  if (i + 1 >= t_.size()) i = t_.size() - 2;
  auto jt = std::upper_bound(p_.begin(), p_.end(), p);
  std::size_t j = static_cast<std::size_t>(std::distance(p_.begin(), jt));
  j = (j == 0) ? 0 : j - 1;
  if (j + 1 >= np) j = np - 2;

  const double ft = (T - t_[i]) / (t_[i + 1] - t_[i]);
  const double fp = (std::log(p) - lnp_[j]) / (lnp_[j + 1] - lnp_[j]);
  const double v00 = field[i * np + j];
  const double v01 = field[i * np + j + 1];
  const double v10 = field[(i + 1) * np + j];
  const double v11 = field[(i + 1) * np + j + 1];
  if (!std::isfinite(v00) || !std::isfinite(v01) || !std::isfinite(v10) || !std::isfinite(v11)) {
    std::ostringstream os;
    os << "coolant " << name_ << ": the property table has no valid data around T = " << T
       << " K, p = " << p << " Pa (the state is inside or too close to the two-phase dome)";
    throw RangeError(os.str());
  }
  return (1 - ft) * ((1 - fp) * v00 + fp * v01) + ft * ((1 - fp) * v10 + fp * v11);
}

double CoolantFluid::saturationTemperature(double p) const {
  if (p >= p_crit_) {
    std::ostringstream os;
    os << "coolant " << name_ << ": no saturation temperature above the critical pressure ("
       << p << " Pa >= " << p_crit_ << " Pa)";
    throw RangeError(os.str());
  }
  if (p < p_.front())
    throw RangeError("coolant " + name_ + ": pressure below the tabulated saturation curve");
  std::size_t j = 0;
  while (j + 2 < p_.size() && p_[j + 1] <= p) ++j;
  if (sat_t_[j] < 0.0 || sat_t_[j + 1] < 0.0) {
    // The upper node is supercritical; fall back to the lower one.
    return sat_t_[j] > 0.0 ? sat_t_[j] : sat_t_[j + 1];
  }
  const double f = (std::log(p) - lnp_[j]) / (lnp_[j + 1] - lnp_[j]);
  return (1 - f) * sat_t_[j] + f * sat_t_[j + 1];
}

void CoolantFluid::validTemperatureRange(double p, double& t_lo, double& t_hi) const {
  std::size_t j = 0;
  if (p > p_.front()) {
    auto it = std::upper_bound(p_.begin(), p_.end(), p);
    j = static_cast<std::size_t>(std::distance(p_.begin(), it));
    j = (j == 0) ? 0 : j - 1;
    if (j + 1 >= p_.size()) j = p_.size() - 2;
  }
  // Interpolation touches both bracketing columns, so take the intersection.
  const int lo = std::max(first_valid_[j], first_valid_[j + 1]);
  const int hi = std::min(last_valid_[j], last_valid_[j + 1]);
  t_lo = t_[static_cast<std::size_t>(lo)];
  t_hi = t_[static_cast<std::size_t>(hi)];
}

CoolantState CoolantFluid::at(double T, double p) const {
  if (p < p_.front() || p > p_.back()) {
    std::ostringstream os;
    os << "coolant " << name_ << ": pressure " << p << " Pa outside the tabulated range ["
       << p_.front() << ", " << p_.back() << "] Pa";
    throw RangeError(os.str());
  }
  double t_lo = 0.0, t_hi = 0.0;
  validTemperatureRange(p, t_lo, t_hi);
  if (T < t_lo || T > t_hi) {
    std::ostringstream os;
    os << "coolant " << name_ << ": temperature " << T << " K is outside the valid range ["
       << t_lo << ", " << t_hi << "] K at " << p * 1e-6 << " MPa";
    if (T < t_lo && t_lo > t_.front())
      os << ". The lower limit rises with pressure because the reference equation of state "
            "has no fluid solution below the melting line -- the coolant would freeze.";
    throw RangeError(os.str());
  }
  CoolantState st;
  st.T = T;
  st.p = p;
  st.rho = interp(rho_, T, p);
  st.cp = interp(cp_, T, p);
  st.h = interp(h_, T, p);
  st.mu = interp(mu_, T, p);
  st.k = interp(k_, T, p);
  st.prandtl = st.mu * st.cp / st.k;
  st.supercritical_pressure = (p >= p_crit_);
  if (st.supercritical_pressure) {
    st.phase = CoolantPhase::kSupercriticalPressure;
  } else {
    st.t_saturation = saturationTemperature(p);
    if (T >= t_crit_) {
      // Above the critical temperature there is no liquid branch at all, so
      // the fluid is a dense gas and cannot boil however low the pressure is.
      st.phase = CoolantPhase::kSupercriticalTemperature;
    } else if (T < st.t_saturation) {
      st.phase = CoolantPhase::kLiquid;
      st.near_saturation = (st.t_saturation - T) < CoolantState::kSaturationMargin * st.t_saturation;
    } else {
      st.phase = CoolantPhase::kVapor;
    }
  }
  return st;
}

double CoolantFluid::temperatureFromEnthalpy(double h, double p, double T_guess) const {
  double lo = 0.0, hi = 0.0;
  validTemperatureRange(p, lo, hi);
  const double h_lo = interp(h_, lo, p), h_hi = interp(h_, hi, p);
  if (h < h_lo || h > h_hi) {
    std::ostringstream os;
    os << "coolant " << name_ << ": enthalpy " << h << " J/kg at " << p
       << " Pa is outside the tabulated range [" << h_lo << ", " << h_hi << "] J/kg";
    throw RangeError(os.str());
  }
  double T = std::min(std::max(T_guess, lo), hi);
  for (int i = 0; i < 200; ++i) {
    const double f = interp(h_, T, p) - h;
    if (std::abs(f) < 1.0e-8 * std::max(1.0, std::abs(h))) return T;
    if (f > 0.0) hi = T; else lo = T;
    const double cp = interp(cp_, T, p);
    double T_new = (cp > 0.0) ? T - f / cp : 0.5 * (lo + hi);
    if (!(T_new > lo && T_new < hi)) T_new = 0.5 * (lo + hi);
    if (std::abs(T_new - T) < 1e-12 * T) return T_new;
    T = T_new;
  }
  throw ConvergenceError("coolant " + name_ + ": temperature-from-enthalpy did not converge");
}

// ---------------------------------------------------------------------------
// Peng-Robinson
// ---------------------------------------------------------------------------

PengRobinsonFluid::PengRobinsonFluid(double t_crit, double p_crit, double acentric,
                                     double molar_mass)
    : t_crit_(t_crit), p_crit_(p_crit), acentric_(acentric), molar_mass_(molar_mass) {
  if (!(t_crit_ > 0.0 && p_crit_ > 0.0 && molar_mass_ > 0.0))
    throw ConfigError("Peng-Robinson: critical constants must be positive");
  const double R = constants::R_universal;
  a_c_ = 0.45724 * R * R * t_crit_ * t_crit_ / p_crit_;
  b_ = 0.07780 * R * t_crit_ / p_crit_;
  kappa_ = 0.37464 + 1.54226 * acentric_ - 0.26992 * acentric_ * acentric_;
}

double PengRobinsonFluid::compressibility(double T, double p) const {
  const double R = constants::R_universal;
  const double u = 1.0 + kappa_ * (1.0 - std::sqrt(T / t_crit_));
  const double a = a_c_ * u * u;
  const double A = a * p / (R * R * T * T);
  const double B = b_ * p / (R * T);
  // Z^3 + c2 Z^2 + c1 Z + c0 = 0
  const double c2 = -(1.0 - B);
  const double c1 = A - 3.0 * B * B - 2.0 * B;
  const double c0 = -(A * B - B * B - B * B * B);
  // Solve by Newton from both ends and keep the physical roots.
  auto f = [&](double Z) { return ((Z + c2) * Z + c1) * Z + c0; };
  auto df = [&](double Z) { return (3.0 * Z + 2.0 * c2) * Z + c1; };
  double best = 1.0;
  bool found = false;
  for (double Z0 : {B * 1.0001 + 1e-9, 0.3, 1.0, 3.0}) {
    double Z = Z0;
    bool ok = false;
    for (int i = 0; i < 200; ++i) {
      const double d = df(Z);
      if (d == 0.0) break;
      const double Zn = Z - f(Z) / d;
      if (!std::isfinite(Zn)) break;
      if (std::abs(Zn - Z) < 1e-14) { Z = Zn; ok = true; break; }
      Z = Zn;
    }
    if (!ok || Z <= B) continue;
    // Prefer the root with the lower molar Gibbs energy; for a single phase
    // there is only one real root anyway.
    if (!found || Z > best) { best = Z; found = true; }
  }
  if (!found) throw ConvergenceError("Peng-Robinson: no physical compressibility root");
  return best;
}

double PengRobinsonFluid::density(double T, double p) const {
  const double Z = compressibility(T, p);
  return p * molar_mass_ / (Z * constants::R_universal * T);
}

double PengRobinsonFluid::enthalpyDeparture(double T, double p) const {
  const double R = constants::R_universal;
  const double u = 1.0 + kappa_ * (1.0 - std::sqrt(T / t_crit_));
  const double a = a_c_ * u * u;
  const double dadT = -a_c_ * kappa_ * u / std::sqrt(T * t_crit_);
  const double Z = compressibility(T, p);
  const double B = b_ * p / (R * T);
  const double s2 = std::sqrt(2.0);
  const double lg = std::log((Z + (1.0 + s2) * B) / (Z + (1.0 - s2) * B));
  const double dh_molar = R * T * (Z - 1.0) + (T * dadT - a) / (2.0 * s2 * b_) * lg;
  return dh_molar / molar_mass_;
}

double PengRobinsonFluid::cpDeparture(double T, double p) const {
  const double R = constants::R_universal;
  const double u = 1.0 + kappa_ * (1.0 - std::sqrt(T / t_crit_));
  const double a = a_c_ * u * u;
  const double dadT = -a_c_ * kappa_ * u / std::sqrt(T * t_crit_);
  const double d2adT2 = a_c_ * kappa_ * (1.0 + kappa_) / (2.0 * std::sqrt(T * T * T * t_crit_));
  const double Z = compressibility(T, p);
  const double v = Z * R * T / p;                 // molar volume
  const double s2 = std::sqrt(2.0);
  const double lg = std::log((v + (1.0 - s2) * b_) / (v + (1.0 + s2) * b_));
  const double dcv = d2adT2 / (2.0 * s2 * b_) * lg;                 // cv - cv_ideal
  const double den = v * v + 2.0 * b_ * v - b_ * b_;
  const double dpdT = R / (v - b_) - dadT / den;
  const double dpdv = -R * T / ((v - b_) * (v - b_)) + a * (2.0 * v + 2.0 * b_) / (den * den);
  const double cp_minus_cv = -T * dpdT * dpdT / dpdv;               // molar
  const double dcp_molar = dcv + cp_minus_cv - R;
  return dcp_molar / molar_mass_;
}

}  // namespace ignis
