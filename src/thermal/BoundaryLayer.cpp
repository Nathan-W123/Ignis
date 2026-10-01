// SPDX-License-Identifier: MIT
#include "ignis/thermal/BoundaryLayer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

#include "ignis/core/Exceptions.hpp"

namespace ignis {

double eckertReferenceTemperature(double t_edge, double t_wall, double t_adiabatic_wall) {
  return t_edge + 0.5 * (t_wall - t_edge) + 0.22 * (t_adiabatic_wall - t_edge);
}

double compressibleShapeFactor(double tw_te, double taw_te, double exponent) {
  if (!(tw_te > 0.0) || !(taw_te > 0.0) || !(exponent > 0.0))
    throw ConfigError("boundary layer: temperature ratios and the profile exponent must be positive");
  // The integrals are taken in u = (y/delta)^(1/n) rather than in y/delta:
  // dy/delta = n u^(n-1) du, which removes the y^(1/n) cusp at the wall and
  // leaves a smooth integrand (a polynomial over the quadratic Crocco-Busemann
  // temperature).  The midpoint rule is then second order; 400 points put H
  // within 2e-5 of the exact (n+2)/n in the incompressible limit.
  constexpr int kPoints = 400;
  double dstar = 0.0, theta = 0.0;
  for (int k = 0; k < kPoints; ++k) {
    const double u = (k + 0.5) / kPoints;
    const double weight = exponent * std::pow(u, exponent - 1.0);
    const double t = tw_te + (taw_te - tw_te) * u - (taw_te - 1.0) * u * u;
    const double rho = 1.0 / t;
    dstar += (1.0 - rho * u) * weight;
    theta += rho * u * (1.0 - u) * weight;
  }
  if (!(theta > 0.0)) throw ConvergenceError("boundary layer: non-positive momentum integral");
  return dstar / theta;
}

double boundaryLayerFilmCoefficient(const BoundaryLayerEdge& e, double enthalpy_thickness,
                                    double t_wall, const ReferenceState& ref) {
  (void)t_wall;  // enters through the reference state the caller evaluated at it
  if (!(enthalpy_thickness > 0.0))
    throw ConfigError("boundary layer: the enthalpy thickness must be positive");
  const double G = e.rho * e.u;
  const double re = G * enthalpy_thickness / ref.viscosity;
  const double st = 0.0128 * (ref.rho / e.rho) * std::pow(re, -0.25) *
                    std::pow(ref.prandtl, -2.0 / 3.0);
  return st * G * ref.cp;
}

BoundaryLayerPoint BoundaryLayerSolution::at(double x) const {
  if (points.empty()) throw ConfigError("boundary layer: empty solution");
  if (x <= points.front().x) return points.front();
  if (x >= points.back().x) return points.back();
  const auto it = std::upper_bound(points.begin(), points.end(), x,
                                   [](double v, const BoundaryLayerPoint& p) { return v < p.x; });
  const auto& b = *it;
  const auto& a = *(it - 1);
  const double t = (x - a.x) / (b.x - a.x);
  auto lerp = [t](double va, double vb) { return va + t * (vb - va); };
  BoundaryLayerPoint p;
  p.x = x;
  p.momentum_thickness = lerp(a.momentum_thickness, b.momentum_thickness);
  p.enthalpy_thickness = lerp(a.enthalpy_thickness, b.enthalpy_thickness);
  p.displacement_thickness = lerp(a.displacement_thickness, b.displacement_thickness);
  p.shape_factor = lerp(a.shape_factor, b.shape_factor);
  p.skin_friction = lerp(a.skin_friction, b.skin_friction);
  p.stanton = lerp(a.stanton, b.stanton);
  p.h_gas = lerp(a.h_gas, b.h_gas);
  p.t_reference = lerp(a.t_reference, b.t_reference);
  p.re_theta = lerp(a.re_theta, b.re_theta);
  p.acceleration = lerp(a.acceleration, b.acceleration);
  return p;
}

namespace {

/// Derivative of a sampled function by second-order differences on a
/// non-uniform grid (one-sided at the ends).
std::vector<double> derivative(const std::vector<double>& x, const std::vector<double>& f) {
  const std::size_t n = x.size();
  std::vector<double> d(n, 0.0);
  if (n < 2) return d;
  d[0] = (f[1] - f[0]) / (x[1] - x[0]);
  d[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);
  for (std::size_t i = 1; i + 1 < n; ++i) {
    const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
    d[i] = (-h1 / (h0 * (h0 + h1))) * f[i - 1] + ((h1 - h0) / (h0 * h1)) * f[i] +
           (h0 / (h1 * (h0 + h1))) * f[i + 1];
  }
  return d;
}

}  // namespace

BoundaryLayerSolution marchBoundaryLayer(const std::vector<BoundaryLayerEdge>& edge,
                                         const ReferenceProperties& properties,
                                         const BoundaryLayerOptions& opts) {
  const std::size_t n = edge.size();
  if (n < 3) throw ConfigError("boundary layer: at least three edge stations are needed");
  for (std::size_t i = 0; i < n; ++i) {
    const auto& e = edge[i];
    if (i > 0 && !(e.x > edge[i - 1].x))
      throw ConfigError("boundary layer: edge stations must be strictly increasing in x");
    if (!(e.radius > 0.0 && e.u > 0.0 && e.rho > 0.0 && e.T > 0.0 && e.p > 0.0))
      throw ConfigError("boundary layer: every edge station needs positive r, u, rho, T and p");
    if (!(e.t_adiabatic_wall > e.t_wall))
      throw ConfigError("boundary layer: the wall must be cooler than the adiabatic wall "
                        "temperature for the energy integral to be defined");
  }
  if (opts.initial_momentum_thickness < 0.0 || opts.initial_enthalpy_thickness < 0.0)
    throw ConfigError("boundary layer: initial thicknesses cannot be negative");
  if (!(opts.start_re_theta > 0.0))
    throw ConfigError("boundary layer: the starting Reynolds number must be positive");
  if (opts.substeps < 1) throw ConfigError("boundary layer: substeps must be at least 1");

  // Reference properties and the driving enthalpy difference at every station.
  std::vector<ReferenceState> ref(n);
  std::vector<double> x(n), ln_u(n), ln_r(n), ln_gdh(n), shape(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& e = edge[i];
    const double ts = eckertReferenceTemperature(e.T, e.t_wall, e.t_adiabatic_wall);
    ref[i] = properties(i, ts);
    x[i] = e.x;
    ln_u[i] = std::log(e.u);
    ln_r[i] = std::log(e.radius);
    ln_gdh[i] = std::log(e.rho * e.u * ref[i].cp * (e.t_adiabatic_wall - e.t_wall));
    shape[i] = compressibleShapeFactor(e.t_wall / e.T, e.t_adiabatic_wall / e.T,
                                       opts.profile_exponent);
  }
  const auto d_ln_u = derivative(x, ln_u);
  const auto d_ln_r = derivative(x, ln_r);
  const auto d_ln_gdh = derivative(x, ln_gdh);

  // Within an interval every coefficient is interpolated linearly; the
  // closure's thickness dependence is evaluated exactly at each RK stage.
  struct Coeffs {
    double G, rho_e, mach, H, a_mom, a_en, rho_ref, mu_ref, pr_ref;
  };
  auto coeffs = [&](std::size_t i, double t) {
    auto lerp = [t](double a, double b) { return a + t * (b - a); };
    const auto& a = edge[i];
    const auto& b = edge[i + 1];
    Coeffs c{};
    c.G = lerp(a.rho * a.u, b.rho * b.u);
    c.rho_e = lerp(a.rho, b.rho);
    c.mach = lerp(a.mach, b.mach);
    c.H = lerp(shape[i], shape[i + 1]);
    const double dlu = lerp(d_ln_u[i], d_ln_u[i + 1]);
    const double dlr = lerp(d_ln_r[i], d_ln_r[i + 1]);
    const double dlg = lerp(d_ln_gdh[i], d_ln_gdh[i + 1]);
    c.a_mom = (2.0 + c.H - c.mach * c.mach) * dlu + dlr;
    c.a_en = dlg + dlr;
    c.rho_ref = lerp(ref[i].rho, ref[i + 1].rho);
    c.mu_ref = lerp(ref[i].viscosity, ref[i + 1].viscosity);
    c.pr_ref = lerp(ref[i].prandtl, ref[i + 1].prandtl);
    return c;
  };
  auto rhs = [](const Coeffs& c, const std::array<double, 2>& y) {
    const double theta = std::max(y[0], 1.0e-12);
    const double delta = std::max(y[1], 1.0e-12);
    const double ratio = c.rho_ref / c.rho_e;
    const double cf2 = 0.0128 * ratio * std::pow(c.G * theta / c.mu_ref, -0.25);
    const double st = 0.0128 * ratio * std::pow(c.G * delta / c.mu_ref, -0.25) *
                      std::pow(c.pr_ref, -2.0 / 3.0);
    return std::array<double, 2>{cf2 - theta * c.a_mom, st - delta * c.a_en};
  };

  BoundaryLayerSolution sol;
  sol.points.resize(n);
  // The thinnest turbulent layer at the first station, unless told otherwise.
  const double theta_min = opts.start_re_theta * ref[0].viscosity / (edge[0].rho * edge[0].u);
  const double theta0 =
      opts.initial_momentum_thickness > 0.0 ? opts.initial_momentum_thickness : theta_min;
  const double delta0 =
      opts.initial_enthalpy_thickness > 0.0 ? opts.initial_enthalpy_thickness : theta0;
  std::array<double, 2> y{theta0, delta0};
  for (std::size_t i = 0; i < n; ++i) {
    if (i > 0) {
      const double hstep = (x[i] - x[i - 1]) / opts.substeps;
      for (int s = 0; s < opts.substeps; ++s) {
        const double t0 = static_cast<double>(s) / opts.substeps;
        const double th = (s + 0.5) / opts.substeps;
        const double t1 = static_cast<double>(s + 1) / opts.substeps;
        const Coeffs c0 = coeffs(i - 1, t0), ch = coeffs(i - 1, th), c1 = coeffs(i - 1, t1);
        const auto k1 = rhs(c0, y);
        const auto k2 = rhs(ch, {y[0] + 0.5 * hstep * k1[0], y[1] + 0.5 * hstep * k1[1]});
        const auto k3 = rhs(ch, {y[0] + 0.5 * hstep * k2[0], y[1] + 0.5 * hstep * k2[1]});
        const auto k4 = rhs(c1, {y[0] + hstep * k3[0], y[1] + hstep * k3[1]});
        for (int j = 0; j < 2; ++j)
          y[j] = std::max(y[j] + hstep / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]),
                          1.0e-12);
      }
    }
    const auto& e = edge[i];
    BoundaryLayerPoint& pt = sol.points[i];
    pt.x = e.x;
    pt.momentum_thickness = y[0];
    pt.enthalpy_thickness = y[1];
    pt.shape_factor = shape[i];
    pt.displacement_thickness = shape[i] * y[0];
    const double G = e.rho * e.u;
    const double ratio = ref[i].rho / e.rho;
    pt.skin_friction = 2.0 * 0.0128 * ratio * std::pow(G * y[0] / ref[i].viscosity, -0.25);
    pt.stanton = 0.0128 * ratio * std::pow(G * y[1] / ref[i].viscosity, -0.25) *
                 std::pow(ref[i].prandtl, -2.0 / 3.0);
    pt.h_gas = pt.stanton * G * ref[i].cp;
    pt.t_reference = eckertReferenceTemperature(e.T, e.t_wall, e.t_adiabatic_wall);
    pt.re_theta = G * y[0] / ref[i].viscosity;
    if (e.viscosity > 0.0)
      pt.acceleration = e.viscosity / (e.rho * e.u * e.u) * d_ln_u[i] * e.u;
  }

  // Validity of the turbulent closure.
  sol.min_re_theta = 1.0e300;
  for (const auto& pt : sol.points) {
    if (pt.acceleration > sol.max_acceleration) {
      sol.max_acceleration = pt.acceleration;
      sol.max_acceleration_x = pt.x;
    }
    sol.min_re_theta = std::min(sol.min_re_theta, pt.re_theta);
  }
  if (sol.max_acceleration > opts.laminarisation_acceleration) {
    std::ostringstream os;
    os << "boundary layer: acceleration parameter K reaches " << sol.max_acceleration
       << " at x = " << sol.max_acceleration_x * 1e3 << " mm, above the " << opts.laminarisation_acceleration
       << " at which a turbulent layer starts to laminarise; the turbulent closure "
          "over-predicts heat transfer there";
    sol.warnings.push_back(os.str());
  }
  // The default start sits exactly on the minimum; round-off must not flag it.
  if (sol.min_re_theta < opts.minimum_turbulent_re_theta * (1.0 - 1.0e-9)) {
    std::ostringstream os;
    os << "boundary layer: momentum-thickness Reynolds number falls to " << sol.min_re_theta
       << ", below the " << opts.minimum_turbulent_re_theta
       << " at which a turbulent closure is meaningful";
    sol.warnings.push_back(os.str());
  }
  return sol;
}

}  // namespace ignis
