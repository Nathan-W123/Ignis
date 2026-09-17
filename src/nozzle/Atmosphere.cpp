// SPDX-License-Identifier: MIT
#include "ignis/nozzle/Atmosphere.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {
namespace {

// Constants *as defined by* the U.S. Standard Atmosphere, 1976.
constexpr double kRstar = 8.31432;      // J/(mol K)  -- the 1976 value
constexpr double kM0 = 28.9644e-3;      // kg/mol
constexpr double kG0 = 9.80665;         // m/s^2
constexpr double kR0 = 6356766.0;       // m, effective Earth radius
constexpr double kGamma = 1.40;         // ratio of specific heats for air
constexpr double kRair = kRstar / kM0;  // J/(kg K)

struct Layer {
  double H_b;   // geopotential base altitude, m'
  double L_b;   // molecular-scale temperature gradient, K/m'
  double T_b;   // base temperature, K
};

// Table 4 of the standard.
constexpr std::array<Layer, 7> kLayers{{
    {     0.0, -6.5e-3, 288.15},
    { 11000.0,  0.0e-3, 216.65},
    { 20000.0,  1.0e-3, 216.65},
    { 32000.0,  2.8e-3, 228.65},
    { 47000.0,  0.0e-3, 270.65},
    { 51000.0, -2.8e-3, 270.65},
    { 71000.0, -2.0e-3, 214.65},
}};
constexpr double kHTop = 84852.0;  // m', top of the 7th layer

/// Geopotential altitude from geometric altitude, m' from m.  The standard's
/// layer table is indexed by geopotential altitude while every caller works in
/// geometric altitude, so the conversion happens once, here.
double geopotential(double z) { return kR0 * z / (kR0 + z); }

/// Base pressures, integrated upward from 101325 Pa at sea level.
const std::array<double, 8>& basePressures() {
  static const std::array<double, 8> p = [] {
    std::array<double, 8> out{};
    out[0] = 101325.0;
    for (std::size_t i = 0; i < kLayers.size(); ++i) {
      const double H_top = (i + 1 < kLayers.size()) ? kLayers[i + 1].H_b : kHTop;
      const double dH = H_top - kLayers[i].H_b;
      if (kLayers[i].L_b == 0.0) {
        out[i + 1] = out[i] * std::exp(-kG0 * kM0 * dH / (kRstar * kLayers[i].T_b));
      } else {
        const double T_top = kLayers[i].T_b + kLayers[i].L_b * dH;
        out[i + 1] = out[i] * std::pow(kLayers[i].T_b / T_top,
                                       kG0 * kM0 / (kRstar * kLayers[i].L_b));
      }
    }
    return out;
  }();
  return p;
}

}  // namespace

AtmosphereState Atmosphere::at(double altitude_m) {
  if (!(altitude_m > -5001.0)) {
    std::ostringstream os;
    os << "atmosphere: altitude " << altitude_m << " m is below the -5000 m model floor";
    throw RangeError(os.str());
  }
  AtmosphereState st;
  st.altitude = altitude_m;

  const double H = geopotential(altitude_m);
  const auto& pb = basePressures();

  if (H <= kHTop) {
    std::size_t i = 0;
    for (std::size_t k = 0; k < kLayers.size(); ++k)
      if (H >= kLayers[k].H_b) i = k;
    const double dH = H - kLayers[i].H_b;
    const double T = kLayers[i].T_b + kLayers[i].L_b * dH;
    double p;
    if (kLayers[i].L_b == 0.0) {
      p = pb[i] * std::exp(-kG0 * kM0 * dH / (kRstar * kLayers[i].T_b));
    } else {
      p = pb[i] * std::pow(kLayers[i].T_b / T, kG0 * kM0 / (kRstar * kLayers[i].L_b));
    }
    st.temperature = T;
    st.pressure = p;
  } else {
    // Exponential extrapolation above the top of the 1976 layer table.
    const double T_top = kLayers.back().T_b + kLayers.back().L_b * (kHTop - kLayers.back().H_b);
    const double scale = kRstar * T_top / (kG0 * kM0);
    st.temperature = T_top;
    st.pressure = pb[7] * std::exp(-(H - kHTop) / scale);
    st.extrapolated = true;
  }
  st.density = st.pressure / (kRair * st.temperature);
  st.speed_of_sound = std::sqrt(kGamma * kRair * st.temperature);
  return st;
}

double Atmosphere::altitudeForPressure(double p) {
  if (!(p > 0.0)) throw RangeError("atmosphere: pressure must be positive");
  const double p_floor = at(-5000.0).pressure;
  if (p > p_floor) {
    std::ostringstream os;
    os << "atmosphere: pressure " << p << " Pa exceeds the -5000 m value of " << p_floor << " Pa";
    throw RangeError(os.str());
  }
  double lo = -5000.0, hi = 1.0e6;
  for (int i = 0; i < 300; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (at(mid).pressure > p) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace ignis
