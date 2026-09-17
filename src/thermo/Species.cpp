// SPDX-License-Identifier: MIT
#include "ignis/thermo/Species.hpp"

#include <cmath>
#include <sstream>
#include <utility>

namespace ignis {

Species::Species(std::string name, std::map<std::string, int> composition, double molar_mass,
                 double t_min, double t_mid, double t_max, std::array<double, 7> low,
                 std::array<double, 7> high, TransportData transport, std::string source)
    : name_(std::move(name)),
      composition_(std::move(composition)),
      molar_mass_(molar_mass),
      t_min_(t_min),
      t_mid_(t_mid),
      t_max_(t_max),
      low_(low),
      high_(high),
      transport_(transport),
      source_(std::move(source)) {
  if (molar_mass_ <= 0.0) throw ConfigError("species " + name_ + ": non-positive molar mass");
  if (!(t_min_ < t_mid_ && t_mid_ < t_max_))
    throw ConfigError("species " + name_ + ": temperature ranges must satisfy Tmin < Tmid < Tmax");
}

int Species::atoms(const std::string& element) const {
  auto it = composition_.find(element);
  return it == composition_.end() ? 0 : it->second;
}

void Species::requireInRange(double T) const {
  if (!inRange(T)) {
    std::ostringstream os;
    os << "species " << name_ << ": temperature " << T << " K outside polynomial range ["
       << t_min_ << ", " << t_max_ << "] K";
    throw RangeError(os.str());
  }
}

void Species::reduced(double T, double& cp_R, double& h_RT, double& s_R) const {
  const auto& a = coeffs(T);
  const double T2 = T * T, T3 = T2 * T, T4 = T3 * T;
  cp_R = a[0] + a[1] * T + a[2] * T2 + a[3] * T3 + a[4] * T4;
  h_RT = a[0] + a[1] * T / 2.0 + a[2] * T2 / 3.0 + a[3] * T3 / 4.0 + a[4] * T4 / 5.0 + a[5] / T;
  s_R = a[0] * std::log(T) + a[1] * T + a[2] * T2 / 2.0 + a[3] * T3 / 3.0 + a[4] * T4 / 4.0 + a[6];
}

double Species::cp(double T) const {
  const auto& a = coeffs(T);
  const double T2 = T * T;
  return constants::R_universal *
         (a[0] + a[1] * T + a[2] * T2 + a[3] * T2 * T + a[4] * T2 * T2);
}

double Species::h(double T) const {
  const auto& a = coeffs(T);
  const double T2 = T * T, T3 = T2 * T, T4 = T3 * T;
  return constants::R_universal * T *
         (a[0] + a[1] * T / 2.0 + a[2] * T2 / 3.0 + a[3] * T3 / 4.0 + a[4] * T4 / 5.0 + a[5] / T);
}

double Species::s0(double T) const {
  const auto& a = coeffs(T);
  const double T2 = T * T, T3 = T2 * T, T4 = T3 * T;
  return constants::R_universal *
         (a[0] * std::log(T) + a[1] * T + a[2] * T2 / 2.0 + a[3] * T3 / 3.0 + a[4] * T4 / 4.0 + a[6]);
}

double Species::mu(double T, double X, double p) const {
  if (!(X > 0.0)) throw RangeError("species " + name_ + ": chemical potential needs X > 0");
  return g0(T) + constants::R_universal * T * std::log(X * p / constants::p_reference);
}

SpeciesState Species::state(double T) const {
  double cp_R, h_RT, s_R;
  reduced(T, cp_R, h_RT, s_R);
  SpeciesState st;
  st.cp = constants::R_universal * cp_R;
  st.h = constants::R_universal * T * h_RT;
  st.s0 = constants::R_universal * s_R;
  st.g0 = st.h - T * st.s0;
  return st;
}

}  // namespace ignis
