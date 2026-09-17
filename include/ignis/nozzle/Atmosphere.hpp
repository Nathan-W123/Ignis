// SPDX-License-Identifier: MIT
#pragma once
/// \file Atmosphere.hpp
/// \brief U.S. Standard Atmosphere 1976.
///
/// Reference: "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562,
/// NASA-TM-X-74335, U.S. Government Printing Office, Washington D.C., 1976.
///
/// The seven gradient/isothermal layers up to 84.852 km' geopotential altitude
/// are reproduced exactly, using the constants *defined by that document*
/// (R* = 8.31432 J/(mol K), M0 = 28.9644 kg/kmol, g0 = 9.80665 m/s^2,
/// r0 = 6 356 766 m) rather than modern CODATA values, so that the tabulated
/// pressures are recovered.  Layer base pressures are integrated from the sea
/// level value instead of being transcribed, which removes any possibility of
/// a copying error.
///
/// Above 86 km geometric altitude the model switches to an exponential
/// extrapolation with the scale height of the topmost layer.  Ambient pressure
/// there is below 0.4 Pa, i.e. already vacuum as far as nozzle performance is
/// concerned, and `AtmosphereState::extrapolated` flags it.

#include <string>

namespace ignis {

/// Ambient conditions at one altitude.
struct AtmosphereState {
  double altitude = 0.0;     ///< geometric altitude, m
  double pressure = 0.0;     ///< Pa
  double temperature = 0.0;  ///< K
  double density = 0.0;      ///< kg/m^3
  double speed_of_sound = 0.0;  ///< m/s
  /// True when the state came from the exponential extrapolation above 86 km.
  bool extrapolated = false;
};

/// U.S. Standard Atmosphere 1976.
class Atmosphere {
 public:
  /// Highest geometric altitude represented exactly by the 1976 layers, m.
  static constexpr double kTopOfModel = 86000.0;

  /// Conditions at a geometric altitude (m).  Negative altitudes down to
  /// -5000 m are accepted (the lowest layer is simply extended downwards, as
  /// the standard itself does); below that a RangeError is thrown.
  static AtmosphereState at(double altitude_m);

  /// Geometric altitude at which the ambient pressure equals `p`, m.
  /// Throws RangeError if the pressure is outside the representable range.
  static double altitudeForPressure(double p);
};

}  // namespace ignis
