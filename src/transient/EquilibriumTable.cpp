// SPDX-License-Identifier: MIT
#include "ignis/transient/EquilibriumTable.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <atomic>
#include <mutex>
#include <random>
#include <thread>
#include <sstream>

#include "ignis/nozzle/NozzleFlow.hpp"

namespace ignis {

std::string InterpolationError::summary() const {
  std::ostringstream os;
  os << std::scientific << std::setprecision(3);
  os << "equilibrium table interpolation error over " << samples << " random interior points\n"
     << "  temperature : max " << max_rel_temperature << ", rms " << rms_rel_temperature << "\n"
     << "  pressure    : max " << max_rel_pressure << ", rms " << rms_rel_pressure << "\n"
     << "  molar mass  : max " << max_rel_molar_mass << "\n"
     << "  c*          : max " << max_rel_cstar;
  return os.str();
}

EquilibriumTable EquilibriumTable::build(const EquilibriumSolver& solver,
                                         const PropellantMixture& base, const TableGrid& grid,
                                         bool verbose) {
  if (grid.mr_points < 2 || grid.t_points < 2 || grid.p_points < 2)
    throw ConfigError("equilibrium table: every axis needs at least two points");
  if (!(grid.mr_min > 0.0 && grid.mr_max > grid.mr_min))
    throw ConfigError("equilibrium table: mixture-ratio range is invalid");
  if (!(grid.t_min > 0.0 && grid.t_max > grid.t_min))
    throw ConfigError("equilibrium table: temperature range is invalid");
  if (!(grid.p_min > 0.0 && grid.p_max > grid.p_min))
    throw ConfigError("equilibrium table: pressure range is invalid");

  EquilibriumTable tab;
  tab.grid_ = grid;
  tab.mr_.resize(static_cast<std::size_t>(grid.mr_points));
  tab.t_.resize(static_cast<std::size_t>(grid.t_points));
  tab.p_.resize(static_cast<std::size_t>(grid.p_points));
  tab.lnp_.resize(tab.p_.size());
  for (int i = 0; i < grid.mr_points; ++i)
    tab.mr_[static_cast<std::size_t>(i)] =
        grid.mr_min + (grid.mr_max - grid.mr_min) * i / (grid.mr_points - 1);
  for (int j = 0; j < grid.t_points; ++j)
    tab.t_[static_cast<std::size_t>(j)] =
        grid.t_min + (grid.t_max - grid.t_min) * j / (grid.t_points - 1);
  for (int k = 0; k < grid.p_points; ++k) {
    tab.p_[static_cast<std::size_t>(k)] =
        grid.p_min * std::pow(grid.p_max / grid.p_min,
                              static_cast<double>(k) / (grid.p_points - 1));
    tab.lnp_[static_cast<std::size_t>(k)] = std::log(tab.p_[static_cast<std::size_t>(k)]);
  }

  const std::size_t n3 = tab.mr_.size() * tab.t_.size() * tab.p_.size();
  tab.u_.resize(n3); tab.molar_.resize(n3); tab.gamma_.resize(n3); tab.cv_.resize(n3);
  const std::size_t n2 = tab.mr_.size() * tab.p_.size();
  tab.cstar_.resize(n2); tab.tflame_.resize(n2);
  tab.h_react_.resize(tab.mr_.size());
  tab.q_comb_.resize(tab.mr_.size());

  const auto& db = solver.database();
  const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  const int nthreads = static_cast<int>(std::min<unsigned>(hw, static_cast<unsigned>(grid.mr_points)));
  std::atomic<int> next_row{0};
  std::mutex log_mutex;
  std::vector<std::string> errors;
  auto worker = [&]() {
    for (;;) {
      const int i = next_row.fetch_add(1);
      if (i >= grid.mr_points) return;
      try {
    const auto mix = base.withMixtureRatio(tab.mr_[static_cast<std::size_t>(i)]);
    const Eigen::VectorXd b = mix.elementMoles(db);
    tab.h_react_[static_cast<std::size_t>(i)] = mix.enthalpy(db);

    // Heat of combustion: reactant enthalpy minus the equilibrium product
    // enthalpy, both referred to 298.15 K and 1 bar.
    const auto ref = solver.tp(b, 298.15, constants::p_reference);
    tab.q_comb_[static_cast<std::size_t>(i)] =
        tab.h_react_[static_cast<std::size_t>(i)] - ref.state.h;

    for (int k = 0; k < grid.p_points; ++k) {
      const double p = tab.p_[static_cast<std::size_t>(k)];
      const auto hp = solver.hp(b, tab.h_react_[static_cast<std::size_t>(i)], p);
      const NozzleFlow flow(solver, CompositionModel::kEquilibrium,
                            ChamberReference{b, hp.state});
      const std::size_t i2 = static_cast<std::size_t>(i) * tab.p_.size() +
                             static_cast<std::size_t>(k);
      tab.cstar_[i2] = flow.cStarIdeal();
      tab.tflame_[i2] = hp.state.T;

      for (int j = 0; j < grid.t_points; ++j) {
        const double T = tab.t_[static_cast<std::size_t>(j)];
        const auto st = solver.tp(b, T, p);
        const std::size_t id = tab.idx(i, j, k);
        tab.u_[id] = st.state.u;
        tab.molar_[id] = st.state.M;
        tab.gamma_[id] = st.state.gamma_s;
        tab.cv_[id] = st.state.cv_eff;
      }
    }
      } catch (const std::exception& e) {
        std::lock_guard<std::mutex> lock(log_mutex);
        errors.emplace_back("mixture ratio " +
                            std::to_string(tab.mr_[static_cast<std::size_t>(i)]) + ": " + e.what());
        return;
      }
      if (verbose) {
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << "  equilibrium table: mixture ratio "
                  << tab.mr_[static_cast<std::size_t>(i)] << " done\n";
      }
    }
  };
  {
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(nthreads));
    for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
  }
  if (!errors.empty()) {
    std::ostringstream os;
    os << "equilibrium table: " << errors.size() << " row(s) failed:";
    for (const auto& e : errors) os << "\n  " << e;
    throw ConvergenceError(os.str());
  }
  return tab;
}

namespace {

/// Locate the interval [ax[i], ax[i+1]] containing v.
std::size_t locateInterval(const std::vector<double>& ax, double v) {
  auto it = std::upper_bound(ax.begin(), ax.end(), v);
  std::size_t i = static_cast<std::size_t>(std::distance(ax.begin(), it));
  i = (i == 0) ? 0 : i - 1;
  if (i + 1 >= ax.size()) i = ax.size() - 2;
  return i;
}

/// Cubic Hermite interpolation on a (possibly non-uniform) axis using centred
/// difference slopes, falling back to one-sided slopes at the ends.  This is
/// C1 and fourth-order accurate in the node spacing, against second order for
/// the linear interpolant.
double hermite(const std::vector<double>& ax, std::size_t i, double v,
               const std::array<double, 4>& f) {
  const double h = ax[i + 1] - ax[i];
  const double s = (v - ax[i]) / h;
  const double s2 = s * s, s3 = s2 * s;
  const double h00 = 2 * s3 - 3 * s2 + 1;
  const double h10 = s3 - 2 * s2 + s;
  const double h01 = -2 * s3 + 3 * s2;
  const double h11 = s3 - s2;
  const double m0 = (i == 0) ? (f[2] - f[1]) / h : (f[2] - f[0]) / (ax[i + 1] - ax[i - 1]);
  const double m1 = (i + 2 >= ax.size()) ? (f[2] - f[1]) / h
                                         : (f[3] - f[1]) / (ax[i + 2] - ax[i]);
  return h00 * f[1] + h10 * h * m0 + h01 * f[2] + h11 * h * m1;
}

}  // namespace

double EquilibriumTable::interp3(const std::vector<double>& f, double mr, double T,
                                 double p) const {
  // Tensor-product cubic Hermite in (mixture ratio, temperature) -- the two
  // directions carrying real curvature -- and linear in ln p, where the
  // tabulated functions are close to straight.
  const std::size_t i = locateInterval(mr_, mr);
  const std::size_t j = locateInterval(t_, T);
  const std::size_t k = locateInterval(p_, p);
  const double fp = (std::log(p) - lnp_[k]) / (lnp_[k + 1] - lnp_[k]);

  auto clampIdx = [](std::size_t base, int off, std::size_t n) {
    long v = static_cast<long>(base) + off - 1;
    if (v < 0) v = 0;
    if (v > static_cast<long>(n) - 1) v = static_cast<long>(n) - 1;
    return static_cast<int>(v);
  };

  double out = 0.0;
  for (int dk = 0; dk < 2; ++dk) {
    std::array<double, 4> along_mr{};
    for (int a = 0; a < 4; ++a) {
      const int ii = clampIdx(i, a, mr_.size());
      std::array<double, 4> along_t{};
      for (int b = 0; b < 4; ++b) {
        const int jj = clampIdx(j, b, t_.size());
        along_t[static_cast<std::size_t>(b)] = f[idx(ii, jj, static_cast<int>(k) + dk)];
      }
      along_mr[static_cast<std::size_t>(a)] = hermite(t_, j, T, along_t);
    }
    const double v = hermite(mr_, i, mr, along_mr);
    out += (dk ? fp : 1.0 - fp) * v;
  }
  return out;
}

double EquilibriumTable::interp2(const std::vector<double>& f, double mr, double p) const {
  const std::size_t i = locateInterval(mr_, mr);
  const std::size_t k = locateInterval(p_, p);
  const double fp = (std::log(p) - lnp_[k]) / (lnp_[k + 1] - lnp_[k]);
  const std::size_t np = p_.size();
  auto clampIdx = [&](int off) {
    long v = static_cast<long>(i) + off - 1;
    if (v < 0) v = 0;
    if (v > static_cast<long>(mr_.size()) - 1) v = static_cast<long>(mr_.size()) - 1;
    return static_cast<std::size_t>(v);
  };
  std::array<double, 4> along_mr{};
  for (int a = 0; a < 4; ++a) {
    const std::size_t ii = clampIdx(a);
    along_mr[static_cast<std::size_t>(a)] =
        (1 - fp) * f[ii * np + k] + fp * f[ii * np + k + 1];
  }
  return hermite(mr_, i, mr, along_mr);
}

double EquilibriumTable::interp1mr(const std::vector<double>& f, double mr) const {
  const std::size_t i = locateInterval(mr_, mr);
  auto clampIdx = [&](int off) {
    long v = static_cast<long>(i) + off - 1;
    if (v < 0) v = 0;
    if (v > static_cast<long>(mr_.size()) - 1) v = static_cast<long>(mr_.size()) - 1;
    return static_cast<std::size_t>(v);
  };
  std::array<double, 4> vals{};
  for (int a = 0; a < 4; ++a) vals[static_cast<std::size_t>(a)] = f[clampIdx(a)];
  return hermite(mr_, i, mr, vals);
}

namespace {
void requireInside(double v, double lo, double hi, const char* what) {
  if (v < lo || v > hi) {
    std::ostringstream os;
    os << "equilibrium table: " << what << " = " << v << " outside the tabulated range [" << lo
       << ", " << hi << "]";
    throw RangeError(os.str());
  }
}
}  // namespace

double EquilibriumTable::internalEnergy(double mr, double T, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(T, grid_.t_min, grid_.t_max, "temperature");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp3(u_, mr, T, p);
}
double EquilibriumTable::molarMass(double mr, double T, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(T, grid_.t_min, grid_.t_max, "temperature");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp3(molar_, mr, T, p);
}
double EquilibriumTable::gammaS(double mr, double T, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(T, grid_.t_min, grid_.t_max, "temperature");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp3(gamma_, mr, T, p);
}
double EquilibriumTable::cv(double mr, double T, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(T, grid_.t_min, grid_.t_max, "temperature");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp3(cv_, mr, T, p);
}
double EquilibriumTable::cStar(double mr, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp2(cstar_, mr, p);
}
double EquilibriumTable::flameTemperature(double mr, double p) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  requireInside(p, grid_.p_min, grid_.p_max, "pressure");
  return interp2(tflame_, mr, p);
}
double EquilibriumTable::reactantEnthalpy(double mr) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  return interp1mr(h_react_, mr);
}
double EquilibriumTable::heatOfCombustion(double mr) const {
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");
  return interp1mr(q_comb_, mr);
}

void EquilibriumTable::solveState(double mr, double u, double rho, double& T, double& p) const {
  if (!(rho > 0.0)) throw RangeError("equilibrium table: density must be positive");
  requireInside(mr, grid_.mr_min, grid_.mr_max, "mixture ratio");

  // Bracket on temperature: u is monotonically increasing in T at fixed rho.
  auto pressureAt = [&](double Tq) {
    double pq = rho * constants::R_universal * Tq / 0.022;   // M ~ 22 g/mol first guess
    for (int i = 0; i < 100; ++i) {
      pq = std::min(std::max(pq, grid_.p_min), grid_.p_max);
      const double M = interp3(molar_, mr, Tq, pq);
      const double pn = rho * constants::R_universal * Tq / M;
      const double rel = std::abs(pn - pq) / pn;
      pq = std::min(std::max(pn, grid_.p_min), grid_.p_max);
      if (rel < 1e-13) break;
    }
    return pq;
  };
  auto residual = [&](double Tq) {
    const double pq = pressureAt(Tq);
    return interp3(u_, mr, Tq, pq) - u;
  };

  double lo = grid_.t_min, hi = grid_.t_max;
  double f_lo = residual(lo), f_hi = residual(hi);
  // A state exactly on an axis endpoint must not be rejected by round-off, so
  // the bracket test carries a relative tolerance.
  const double tol = 1.0e-9 * std::max(1.0, std::abs(u));
  if (f_lo > tol || f_hi < -tol) {
    std::ostringstream os;
    os << "equilibrium table: specific internal energy " << u << " J/kg at rho = " << rho
       << " kg/m^3 and O/F = " << mr << " lies outside the tabulated span ["
       << (u + f_lo) << ", " << (u + f_hi) << "] J/kg (temperature axis "
       << grid_.t_min << " to " << grid_.t_max << " K).";
    if (f_lo > 0.0) {
      os << "\nThe chamber energy is below what equilibrium products can hold at "
         << grid_.t_min << " K. In a transient this normally means the run is spending a "
            "long interval with a large propellant lead and almost no heat release: the "
            "zero-dimensional products model represents the chamber contents as combustion "
            "products at all times, so it is valid from ignition onwards, not through a "
            "cold pre-ignition fill. Start the ignition ramp as soon as both propellants "
            "are present, shorten the valve lead, or lower the table temperature floor "
            "(the NASA polynomials stop at 200 K).";
    }
    throw RangeError(os.str());
  }
  if (f_lo >= 0.0) { T = lo; p = pressureAt(lo); return; }
  if (f_hi <= 0.0) { T = hi; p = pressureAt(hi); return; }
  T = 0.5 * (lo + hi);
  for (int i = 0; i < 200; ++i) {
    T = 0.5 * (lo + hi);
    const double f = residual(T);
    if (f > 0.0) hi = T; else lo = T;
    if (hi - lo < 1e-9 * std::max(1.0, T)) break;
  }
  p = pressureAt(T);
}

InterpolationError EquilibriumTable::measureError(const EquilibriumSolver& solver,
                                                  const PropellantMixture& base, int samples,
                                                  unsigned seed) const {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  InterpolationError err;
  err.samples = samples;
  double sum_t = 0.0, sum_p = 0.0;
  const auto& db = solver.database();
  int used = 0;
  // Sample strictly inside the grid: the round trip through solveState is only
  // defined where the bracket has room on both sides, and the quantity of
  // interest is the interior interpolation error.
  const double dt_edge = 0.5 * (grid_.t_max - grid_.t_min) / (grid_.t_points - 1);
  const double dmr_edge = 0.5 * (grid_.mr_max - grid_.mr_min) / (grid_.mr_points - 1);
  for (int s = 0; s < samples; ++s) {
    const double mr = (grid_.mr_min + dmr_edge) +
                      (grid_.mr_max - grid_.mr_min - 2 * dmr_edge) * uni(rng);
    const double T = (grid_.t_min + dt_edge) +
                     (grid_.t_max - grid_.t_min - 2 * dt_edge) * uni(rng);
    const double p = grid_.p_min * std::pow(grid_.p_max / grid_.p_min,
                                            0.02 + 0.96 * uni(rng));
    const auto mix = base.withMixtureRatio(mr);
    const Eigen::VectorXd b = mix.elementMoles(db);
    EquilibriumResult exact;
    try {
      exact = solver.tp(b, T, p);
    } catch (const IgnisError&) {
      continue;
    }
    // Recover (T, p) from the tabulated (u, rho) and compare.
    double Ti = 0.0, pi = 0.0;
    try {
      solveState(mr, exact.state.u, exact.state.rho, Ti, pi);
    } catch (const IgnisError&) {
      continue;  // the round trip left the grid; not an interior sample
    }
    ++used;
    const double et = std::abs(Ti - T) / T;
    const double ep = std::abs(pi - p) / p;
    err.max_rel_temperature = std::max(err.max_rel_temperature, et);
    err.max_rel_pressure = std::max(err.max_rel_pressure, ep);
    sum_t += et * et;
    sum_p += ep * ep;
    err.max_rel_molar_mass = std::max(
        err.max_rel_molar_mass,
        std::abs(molarMass(mr, T, p) - exact.state.M) / exact.state.M);

    const auto hp = solver.hp(b, mix.enthalpy(db), p);
    const NozzleFlow flow(solver, CompositionModel::kEquilibrium, ChamberReference{b, hp.state});
    err.max_rel_cstar = std::max(err.max_rel_cstar,
                                 std::abs(cStar(mr, p) - flow.cStarIdeal()) / flow.cStarIdeal());
  }
  err.samples = used;
  if (used > 0) {
    err.rms_rel_temperature = std::sqrt(sum_t / used);
    err.rms_rel_pressure = std::sqrt(sum_p / used);
  }
  return err;
}

void EquilibriumTable::exportCsv(const std::string& path) const {
  std::ofstream out(path);
  if (!out) throw ConfigError("cannot write equilibrium table '" + path + "'");
  out << "# Ignis equilibrium property table\n";
  out << "#grid mr_min=" << grid_.mr_min << " mr_max=" << grid_.mr_max
      << " mr_points=" << grid_.mr_points << " t_min=" << grid_.t_min << " t_max=" << grid_.t_max
      << " t_points=" << grid_.t_points << " p_min=" << grid_.p_min << " p_max=" << grid_.p_max
      << " p_points=" << grid_.p_points << "\n";
  out << std::setprecision(12);
  out << "mr,T,p,u,molar_mass,gamma_s,cv\n";
  for (int i = 0; i < grid_.mr_points; ++i)
    for (int j = 0; j < grid_.t_points; ++j)
      for (int k = 0; k < grid_.p_points; ++k) {
        const std::size_t id = idx(i, j, k);
        out << mr_[static_cast<std::size_t>(i)] << ',' << t_[static_cast<std::size_t>(j)] << ','
            << p_[static_cast<std::size_t>(k)] << ',' << u_[id] << ',' << molar_[id] << ','
            << gamma_[id] << ',' << cv_[id] << '\n';
      }
}

}  // namespace ignis
