// SPDX-License-Identifier: MIT
#include "ignis/thermal/FilmCooling.hpp"

#include <algorithm>
#include <cmath>

#include "ignis/core/Exceptions.hpp"

namespace ignis {

double FilmSlot::coolantVelocity() const {
  if (!(density > 0.0 && perimeter > 0.0 && slot_height > 0.0))
    throw ConfigError("film cooling: the slot needs a positive density, length and height");
  return mass_flow / (density * perimeter * slot_height);
}

double FilmSlot::diffusivity() const {
  if (!(density > 0.0 && cp > 0.0 && conductivity > 0.0))
    throw ConfigError("film cooling: the coolant needs a positive density, cp and conductivity");
  return conductivity / (density * cp);
}

double hatchPapellVelocityFunction(double r) {
  if (!(r > 0.0)) throw ConfigError("film cooling: the velocity ratio must be positive");
  if (r >= 1.0) return 1.0 + 0.4 * std::atan(r - 1.0);   // eq. (10), radians
  const double inv = 1.0 / r;
  return std::pow(inv, 1.5 * (inv - 1.0));                // eq. (11)
}

double hatchPapellFactor(const FilmSlot& s, bool* clamped) {
  if (!(s.mass_flow > 0.0)) throw ConfigError("film cooling: the film mass flow must be positive");
  if (!(s.gas_velocity > 0.0)) throw ConfigError("film cooling: the gas velocity must be positive");
  const double vc = s.coolantVelocity();
  const double r = s.gas_velocity / vc;
  const double r_used =
      std::min(std::max(r, kHatchPapellMinVelocityRatio), kHatchPapellMaxVelocityRatio);
  if (clamped != nullptr) *clamped = (r_used != r);
  // Eq. (12) uses the gas velocity in S V/alpha_c; f carries the ratio.
  return std::pow(s.slot_height * s.gas_velocity / s.diffusivity(), 0.125) *
         hatchPapellVelocityFunction(r_used);
}

double hatchPapellEffectiveness(double heat_group, double factor) {
  if (heat_group < 0.0) throw ConfigError("film cooling: the heat-capacity group cannot be negative");
  if (!(factor > 0.0)) throw ConfigError("film cooling: the correlation factor must be positive");
  if (heat_group <= 0.04) return 1.0;
  return std::exp(-(heat_group - 0.04) * factor);
}

double hatchPapellGasCoefficient(double density, double velocity, double diameter,
                                 double viscosity, double conductivity, double prandtl) {
  if (!(density > 0.0 && velocity > 0.0 && diameter > 0.0 && viscosity > 0.0 &&
        conductivity > 0.0 && prandtl > 0.0))
    throw ConfigError("film cooling: gas-side coefficient inputs must be positive");
  const double re = density * velocity * diameter / viscosity;
  return 0.0265 * conductivity / diameter * std::pow(re, 0.8) * std::pow(prandtl, 0.3);
}

}  // namespace ignis
