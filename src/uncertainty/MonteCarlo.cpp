// SPDX-License-Identifier: MIT
#include "ignis/uncertainty/MonteCarlo.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <thread>

#include <Eigen/Dense>

#include "ignis/optimize/Parameters.hpp"

namespace ignis {

std::string toString(DistributionType d) {
  switch (d) {
    case DistributionType::kNormal: return "normal";
    case DistributionType::kLogNormal: return "lognormal";
    case DistributionType::kUniform: return "uniform";
    case DistributionType::kTriangular: return "triangular";
  }
  return "?";
}

DistributionType distributionFromString(const std::string& s) {
  if (s == "normal") return DistributionType::kNormal;
  if (s == "lognormal") return DistributionType::kLogNormal;
  if (s == "uniform") return DistributionType::kUniform;
  if (s == "triangular") return DistributionType::kTriangular;
  throw ConfigError("unknown distribution '" + s +
                    "'; use normal, lognormal, uniform or triangular");
}

namespace {

/// SplitMix64 (Vigna).  Used to derive a per-sample seed that does not depend
/// on the thread count or the scheduling order.
std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ull;
  std::uint64_t z = x;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

double drawOnce(const UncertainInput& in, double nominal, std::mt19937_64& rng) {
  double a = in.a, b = in.b, c = in.c;
  if (in.relative) { a *= nominal; b *= nominal; c *= nominal; }
  switch (in.distribution) {
    case DistributionType::kNormal: {
      std::normal_distribution<double> d(a, in.relative ? in.b * std::abs(nominal) : b);
      return d(rng);
    }
    case DistributionType::kLogNormal: {
      // a is the median, in.b the standard deviation of ln(x).
      std::normal_distribution<double> d(std::log(a), in.b);
      return std::exp(d(rng));
    }
    case DistributionType::kUniform: {
      std::uniform_real_distribution<double> d(a, b);
      return d(rng);
    }
    case DistributionType::kTriangular: {
      std::uniform_real_distribution<double> u(0.0, 1.0);
      const double lo = a, mode = b, hi = c;
      if (!(hi > lo) || mode < lo || mode > hi)
        throw ConfigError("triangular distribution needs low <= mode <= high with high > low");
      const double f = (mode - lo) / (hi - lo);
      const double r = u(rng);
      return (r < f) ? lo + std::sqrt(r * (hi - lo) * (mode - lo))
                     : hi - std::sqrt((1.0 - r) * (hi - lo) * (hi - mode));
    }
  }
  throw ConfigError("unhandled distribution");
}

double draw(const UncertainInput& in, double nominal, std::mt19937_64& rng) {
  for (int attempt = 0; attempt < 1000; ++attempt) {
    const double v = drawOnce(in, nominal, rng);
    if (v >= in.lower_bound && v <= in.upper_bound && std::isfinite(v)) return v;
  }
  std::ostringstream os;
  os << "Monte Carlo: could not draw a value for '" << in.parameter
     << "' inside its bounds [" << in.lower_bound << ", " << in.upper_bound
     << "] after 1000 attempts; the bounds and the distribution are inconsistent";
  throw ConfigError(os.str());
}

/// Ranks with average ties, used by the Spearman coefficient.
std::vector<double> ranks(const std::vector<double>& v) {
  const std::size_t n = v.size();
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
  std::vector<double> r(n, 0.0);
  std::size_t i = 0;
  while (i < n) {
    std::size_t j = i;
    while (j + 1 < n && v[idx[j + 1]] == v[idx[i]]) ++j;
    const double avg = 0.5 * (static_cast<double>(i) + static_cast<double>(j)) + 1.0;
    for (std::size_t k = i; k <= j; ++k) r[idx[k]] = avg;
    i = j + 1;
  }
  return r;
}

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = x.size();
  if (n < 2) return 0.0;
  const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
  const double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
  double sxy = 0.0, sxx = 0.0, syy = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double dx = x[i] - mx, dy = y[i] - my;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return 0.0;
  return sxy / std::sqrt(sxx * syy);
}

}  // namespace

Statistics computeStatistics(std::vector<double> v) {
  Statistics s;
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  s.count = static_cast<int>(v.size());
  if (v.empty()) return s;
  std::sort(v.begin(), v.end());
  s.minimum = v.front();
  s.maximum = v.back();
  s.mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  double m2 = 0.0, m3 = 0.0;
  for (double x : v) {
    const double d = x - s.mean;
    m2 += d * d;
    m3 += d * d * d;
  }
  if (v.size() > 1) s.stddev = std::sqrt(m2 / (v.size() - 1));
  if (s.stddev > 0.0) s.skewness = (m3 / v.size()) / std::pow(m2 / v.size(), 1.5);
  auto pct = [&](double q) {
    const double pos = q * (v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(v.size() - 1, lo + 1);
    const double f = pos - lo;
    return (1 - f) * v[lo] + f * v[hi];
  };
  s.p5 = pct(0.05);
  s.p50 = pct(0.50);
  s.p95 = pct(0.95);
  s.p99 = pct(0.99);
  return s;
}

MonteCarloResult runMonteCarlo(const SteadyEngine& engine, const MonteCarloSpec& spec) {
  if (spec.inputs.empty()) throw ConfigError("Monte Carlo: no dispersed inputs");
  if (spec.outputs.empty()) throw ConfigError("Monte Carlo: no outputs requested");
  if (spec.samples < 2) throw ConfigError("Monte Carlo: at least two samples are required");
  for (const auto& in : spec.inputs) (void)parameterUnits(in.parameter);
  for (const auto& o : spec.outputs) (void)metricUnits(o);

  const std::size_t ni = spec.inputs.size();
  const std::size_t no = spec.outputs.size();
  const std::size_t N = static_cast<std::size_t>(spec.samples);

  std::vector<std::vector<double>> xs(ni, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> ys(no, std::vector<double>(N, std::nan("")));
  std::vector<std::string> status(N, "ok");
  std::vector<std::string> messages(N);
  std::vector<double> wall_flag(N, 0.0), boil_flag(N, 0.0), sep_flag(N, 0.0), shock_flag(N, 0.0);
  std::vector<double> stiff_flag(N, 0.0);
  std::vector<double> energy_residual(N, std::nan(""));

  std::vector<double> nominal(ni, 0.0);
  for (std::size_t i = 0; i < ni; ++i)
    nominal[i] = readParameter(engine.config(), spec.inputs[i].parameter);

  const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  const int nthreads = spec.threads > 0 ? spec.threads : static_cast<int>(hw);
  std::atomic<std::size_t> next{0};

  auto worker = [&]() {
    for (;;) {
      const std::size_t k = next.fetch_add(1);
      if (k >= N) return;
      // Per-sample seed: independent of thread count and scheduling order.
      const std::uint64_t s = splitmix64(spec.seed ^ splitmix64(k + 0x9E3779B97F4A7C15ull));
      std::mt19937_64 rng(s);
      EngineConfig cfg = engine.config();
      cfg.sample_profile = false;
      for (std::size_t i = 0; i < ni; ++i) {
        const double v = draw(spec.inputs[i], nominal[i], rng);
        xs[i][k] = v;
        applyParameter(cfg, spec.inputs[i].parameter, v);
      }
      try {
        const auto res = engine.runWith(cfg);
        for (std::size_t o = 0; o < no; ++o) ys[o][k] = readMetric(res, spec.outputs[o]);
        if (res.has_cooling) {
          wall_flag[k] = res.cooling.wall_limit_exceeded ? 1.0 : 0.0;
          boil_flag[k] = res.cooling.boiling_detected ? 1.0 : 0.0;
          energy_residual[k] = res.cooling.energy_balance_residual;
        }
        sep_flag[k] = res.performance.separation_predicted ? 1.0 : 0.0;
        shock_flag[k] = res.performance.shock_in_nozzle ? 1.0 : 0.0;
        if (res.has_feed) stiff_flag[k] = res.feed.stiffness_ok ? 0.0 : 1.0;
      } catch (const std::exception& e) {
        status[k] = "failed";
        messages[k] = e.what();
      }
    }
  };

  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(nthreads));
  for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker);
  for (auto& th : pool) th.join();

  MonteCarloResult res;
  res.requested = spec.samples;
  for (const auto& in : spec.inputs) res.input_names.push_back(in.parameter);
  res.output_names = spec.outputs;

  for (std::size_t k = 0; k < N; ++k) {
    if (status[k] == "ok") {
      ++res.succeeded;
    } else {
      ++res.failed;
      // Collapse the message to its first line so the histogram of failure
      // modes stays readable.
      std::string key = messages[k].substr(0, messages[k].find('\n'));
      if (key.size() > 160) key = key.substr(0, 160) + "...";
      ++res.failure_modes[key];
    }
  }
  auto countFlag = [&](const char* name, const std::vector<double>& f) {
    int n = 0;
    for (std::size_t k = 0; k < N; ++k)
      if (status[k] == "ok" && f[k] > 0.5) ++n;
    if (n > 0) res.flags[name] = n;
  };
  countFlag("wall temperature limit exceeded", wall_flag);
  countFlag("coolant boiling detected", boil_flag);
  countFlag("flow separation predicted", sep_flag);
  countFlag("internal normal shock", shock_flag);
  countFlag("injector stiffness below minimum", stiff_flag);

  if (spec.record_samples) {
    for (std::size_t i = 0; i < ni; ++i)
      res.samples.addColumn(spec.inputs[i].parameter, xs[i],
                            parameterUnits(spec.inputs[i].parameter));
    for (std::size_t o = 0; o < no; ++o)
      res.samples.addColumn(spec.outputs[o], ys[o], metricUnits(spec.outputs[o]));
    res.samples.addColumn("cooling_energy_residual", energy_residual, "-");
    res.samples.addColumn("wall_limit_exceeded", wall_flag, "-");
    res.samples.addColumn("coolant_boiling", boil_flag, "-");
    res.samples.addColumn("separation_predicted", sep_flag, "-");
    res.samples.addColumn("shock_in_nozzle", shock_flag, "-");
    res.samples.addColumn("status", status);
  }

  for (std::size_t o = 0; o < no; ++o)
    res.statistics[spec.outputs[o]] = computeStatistics(ys[o]);
  res.statistics["cooling_energy_residual"] = computeStatistics(energy_residual);

  // --- sensitivity ------------------------------------------------------
  std::vector<std::size_t> good;
  good.reserve(N);
  for (std::size_t k = 0; k < N; ++k)
    if (status[k] == "ok") good.push_back(k);

  // Local finite-difference elasticities at the nominal point.
  std::vector<std::vector<double>> elasticity(no, std::vector<double>(ni, std::nan("")));
  {
    SteadyEngineResult base;
    bool base_ok = true;
    try {
      base = engine.run();
    } catch (const std::exception&) {
      base_ok = false;
    }
    if (base_ok) {
      for (std::size_t i = 0; i < ni; ++i) {
        const double x0 = nominal[i];
        const double dx = spec.fd_step * (std::abs(x0) > 0.0 ? std::abs(x0) : 1.0);
        EngineConfig cp = engine.config(), cm = engine.config();
        cp.sample_profile = false;
        cm.sample_profile = false;
        applyParameter(cp, spec.inputs[i].parameter, x0 + dx);
        applyParameter(cm, spec.inputs[i].parameter, x0 - dx);
        try {
          const auto rp = engine.runWith(cp);
          const auto rm = engine.runWith(cm);
          for (std::size_t o = 0; o < no; ++o) {
            const double yp = readMetric(rp, spec.outputs[o]);
            const double ym = readMetric(rm, spec.outputs[o]);
            const double y0 = readMetric(base, spec.outputs[o]);
            const double dydx = (yp - ym) / (2.0 * dx);
            elasticity[o][i] = (std::abs(y0) > 0.0) ? dydx * x0 / y0 : dydx;
          }
        } catch (const std::exception&) {
          // Leave NaN: a step that leaves the feasible region has no local
          // derivative, and inventing one would be worse than saying so.
        }
      }
    }
  }

  for (std::size_t o = 0; o < no; ++o) {
    OutputSensitivity s;
    s.output = spec.outputs[o];
    s.src.assign(ni, 0.0);
    s.spearman.assign(ni, 0.0);
    s.elasticity = elasticity[o];

    std::vector<std::size_t> rows;
    for (std::size_t k : good)
      if (std::isfinite(ys[o][k])) rows.push_back(k);
    if (rows.size() > ni + 2) {
      const Eigen::Index m = static_cast<Eigen::Index>(rows.size());
      const Eigen::Index p = static_cast<Eigen::Index>(ni);
      Eigen::MatrixXd A(m, p + 1);
      Eigen::VectorXd y(m);
      for (Eigen::Index r = 0; r < m; ++r) {
        A(r, 0) = 1.0;
        for (Eigen::Index c = 0; c < p; ++c)
          A(r, c + 1) = xs[static_cast<std::size_t>(c)][rows[static_cast<std::size_t>(r)]];
        y(r) = ys[o][rows[static_cast<std::size_t>(r)]];
      }
      const Eigen::VectorXd beta = A.colPivHouseholderQr().solve(y);
      const Eigen::VectorXd fit = A * beta;
      const double ybar = y.mean();
      const double ss_tot = (y.array() - ybar).square().sum();
      const double ss_res = (y - fit).squaredNorm();
      s.r_squared = (ss_tot > 0.0) ? 1.0 - ss_res / ss_tot : 0.0;
      const double sy = std::sqrt(ss_tot / std::max<Eigen::Index>(1, m - 1));
      for (std::size_t i = 0; i < ni; ++i) {
        std::vector<double> xi;
        xi.reserve(rows.size());
        for (std::size_t r : rows) xi.push_back(xs[i][r]);
        const auto st = computeStatistics(xi);
        s.src[i] = (sy > 0.0) ? beta(static_cast<Eigen::Index>(i) + 1) * st.stddev / sy : 0.0;
        std::vector<double> yi;
        yi.reserve(rows.size());
        for (std::size_t r : rows) yi.push_back(ys[o][r]);
        s.spearman[i] = pearson(ranks(xi), ranks(yi));
      }
    }
    res.sensitivity.push_back(std::move(s));
  }
  return res;
}

std::string MonteCarloResult::summary() const {
  std::ostringstream os;
  os << "Monte Carlo: " << succeeded << " of " << requested << " samples succeeded ("
     << failed << " failed)\n";
  if (!failure_modes.empty()) {
    os << "  failure modes:\n";
    for (const auto& kv : failure_modes)
      os << "    " << std::setw(6) << kv.second << "  " << kv.first << "\n";
  }
  if (!flags.empty()) {
    os << "  constraint flags among successful samples:\n";
    for (const auto& kv : flags)
      os << "    " << std::setw(6) << kv.second << "  " << kv.first << "\n";
  }
  os << "\n  " << std::left << std::setw(34) << "output" << std::right << std::setw(14) << "mean"
     << std::setw(14) << "stddev" << std::setw(14) << "p5" << std::setw(14) << "p50"
     << std::setw(14) << "p95" << "\n";
  for (const auto& name : output_names) {
    const auto it = statistics.find(name);
    if (it == statistics.end()) continue;
    const auto& s = it->second;
    os << "  " << std::left << std::setw(34) << name << std::right << std::scientific
       << std::setprecision(5) << std::setw(14) << s.mean << std::setw(14) << s.stddev
       << std::setw(14) << s.p5 << std::setw(14) << s.p50 << std::setw(14) << s.p95 << "\n";
  }
  os << "\n  sensitivity (standardised regression coefficients; |SRC| ranked)\n";
  for (const auto& s : sensitivity) {
    os << "    " << s.output << "  (linear model R^2 = " << std::fixed << std::setprecision(4)
       << s.r_squared << ")\n";
    std::vector<std::size_t> order(input_names.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      return std::abs(s.src[a]) > std::abs(s.src[b]);
    });
    for (std::size_t r = 0; r < order.size(); ++r) {
      const std::size_t i = order[r];
      os << "      " << std::left << std::setw(34) << input_names[i] << std::right << std::fixed
         << std::setprecision(4) << std::setw(10) << s.src[i] << "   (Spearman "
         << std::setw(8) << s.spearman[i] << ", elasticity ";
      if (std::isfinite(s.elasticity[i])) os << std::setw(10) << s.elasticity[i];
      else os << std::setw(10) << "n/a";
      os << ")\n";
    }
  }
  return os.str();
}

Table MonteCarloResult::sensitivityTable() const {
  Table t("sensitivity");
  std::vector<std::string> out, in;
  std::vector<double> src, sp, el, r2;
  for (const auto& s : sensitivity)
    for (std::size_t i = 0; i < input_names.size(); ++i) {
      out.push_back(s.output);
      in.push_back(input_names[i]);
      src.push_back(s.src[i]);
      sp.push_back(s.spearman[i]);
      el.push_back(s.elasticity[i]);
      r2.push_back(s.r_squared);
    }
  t.addColumn("output", out);
  t.addColumn("input", in);
  t.addColumn("src", src, "-");
  t.addColumn("spearman", sp, "-");
  t.addColumn("elasticity", el, "-");
  t.addColumn("r_squared", r2, "-");
  return t;
}

Table MonteCarloResult::statisticsTable() const {
  Table t("statistics");
  std::vector<std::string> names;
  std::vector<double> count, mean, sd, mn, mx, p5, p50, p95, p99, skew;
  for (const auto& kv : statistics) {
    names.push_back(kv.first);
    count.push_back(kv.second.count);
    mean.push_back(kv.second.mean);
    sd.push_back(kv.second.stddev);
    mn.push_back(kv.second.minimum);
    mx.push_back(kv.second.maximum);
    p5.push_back(kv.second.p5);
    p50.push_back(kv.second.p50);
    p95.push_back(kv.second.p95);
    p99.push_back(kv.second.p99);
    skew.push_back(kv.second.skewness);
  }
  t.addColumn("metric", names);
  t.addColumn("count", count, "-");
  t.addColumn("mean", mean, "-");
  t.addColumn("stddev", sd, "-");
  t.addColumn("min", mn, "-");
  t.addColumn("max", mx, "-");
  t.addColumn("p5", p5, "-");
  t.addColumn("p50", p50, "-");
  t.addColumn("p95", p95, "-");
  t.addColumn("p99", p99, "-");
  t.addColumn("skewness", skew, "-");
  return t;
}

MonteCarloSpec parseMonteCarlo(const ConfigNode& root) {
  MonteCarloSpec spec;
  const auto m = root["monte_carlo"];
  m.requireOnly({"samples", "seed", "threads", "inputs", "outputs", "record_samples", "fd_step"});
  spec.samples = m["samples"].integer();
  spec.seed = static_cast<std::uint64_t>(m.optional("seed").integer(20260917));
  spec.threads = m.optional("threads").integer(0);
  spec.record_samples = m.optional("record_samples").boolean(true);
  spec.fd_step = m.optional("fd_step").number(1e-8, 0.1, 1.0e-3);
  const auto ins = m["inputs"];
  if (!ins.raw().IsSequence()) throw ConfigError("monte_carlo: 'inputs' must be a sequence");
  for (std::size_t i = 0; i < ins.raw().size(); ++i) {
    ConfigNode n(ins.raw()[i], ins.path() + "[" + std::to_string(i) + "]", "monte_carlo");
    n.requireOnly({"parameter", "distribution", "mean", "sigma", "median", "sigma_log", "low",
                   "high", "mode", "min", "max", "relative"});
    UncertainInput u;
    u.parameter = n["parameter"].text();
    u.distribution = distributionFromString(n.optional("distribution").text("normal"));
    u.relative = n.optional("relative").boolean(false);
    switch (u.distribution) {
      case DistributionType::kNormal:
        u.a = n["mean"].number();
        u.b = n["sigma"].number();
        break;
      case DistributionType::kLogNormal:
        u.a = n["median"].number();
        u.b = n["sigma_log"].number();
        break;
      case DistributionType::kUniform:
        u.a = n["low"].number();
        u.b = n["high"].number();
        break;
      case DistributionType::kTriangular:
        u.a = n["low"].number();
        u.b = n["mode"].number();
        u.c = n["high"].number();
        break;
    }
    u.lower_bound = n.optional("min").number(-std::numeric_limits<double>::infinity());
    u.upper_bound = n.optional("max").number(std::numeric_limits<double>::infinity());
    spec.inputs.push_back(u);
  }
  spec.outputs = m["outputs"].strings();
  return spec;
}

}  // namespace ignis
