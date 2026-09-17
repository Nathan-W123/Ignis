// SPDX-License-Identifier: MIT
#include "ignis/optimize/Sweep.hpp"

#include <atomic>
#include <cmath>
#include <mutex>
#include <sstream>
#include <thread>

#include "ignis/optimize/Parameters.hpp"

namespace ignis {

SweepAxis SweepAxis::linear(const std::string& parameter, double lo, double hi, int points) {
  if (points < 1) throw ConfigError("sweep axis '" + parameter + "': needs at least one point");
  SweepAxis a;
  a.parameter = parameter;
  a.values.reserve(static_cast<std::size_t>(points));
  for (int i = 0; i < points; ++i)
    a.values.push_back(points == 1 ? lo : lo + (hi - lo) * i / (points - 1));
  return a;
}

SweepAxis SweepAxis::logarithmic(const std::string& parameter, double lo, double hi, int points) {
  if (points < 1) throw ConfigError("sweep axis '" + parameter + "': needs at least one point");
  if (!(lo > 0.0 && hi > 0.0))
    throw ConfigError("sweep axis '" + parameter + "': logarithmic spacing needs positive bounds");
  SweepAxis a;
  a.parameter = parameter;
  for (int i = 0; i < points; ++i)
    a.values.push_back(points == 1 ? lo
                                   : lo * std::pow(hi / lo, static_cast<double>(i) / (points - 1)));
  return a;
}

std::string SweepResult::summary() const {
  std::ostringstream os;
  os << "sweep: " << points << " points, " << failures << " failed";
  if (!failure_messages.empty()) {
    os << "\n  first failures:";
    for (std::size_t i = 0; i < failure_messages.size() && i < 3; ++i)
      os << "\n    " << failure_messages[i];
  }
  return os.str();
}

SweepResult runSweep(const SteadyEngine& engine, const SweepSpec& spec) {
  if (spec.axes.empty()) throw ConfigError("sweep: at least one axis is required");
  if (spec.metrics.empty()) throw ConfigError("sweep: at least one metric is required");
  for (const auto& ax : spec.axes) {
    if (ax.values.empty()) throw ConfigError("sweep axis '" + ax.parameter + "' is empty");
    (void)parameterUnits(ax.parameter);  // validates the name
  }
  for (const auto& m : spec.metrics) (void)metricUnits(m);

  std::size_t total = 1;
  for (const auto& ax : spec.axes) total *= ax.values.size();

  std::vector<std::vector<double>> inputs(spec.axes.size(), std::vector<double>(total, 0.0));
  std::vector<std::vector<double>> outputs(spec.metrics.size(),
                                           std::vector<double>(total, std::nan("")));
  std::vector<std::string> status(total, "ok");
  std::vector<std::string> messages(total);

  // Decode a flat index into the grid, last axis varying fastest.
  auto decode = [&](std::size_t flat, std::vector<double>& values) {
    for (std::size_t a = spec.axes.size(); a-- > 0;) {
      const std::size_t n = spec.axes[a].values.size();
      values[a] = spec.axes[a].values[flat % n];
      flat /= n;
    }
  };

  const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  const int nthreads = spec.threads > 0 ? spec.threads : static_cast<int>(hw);
  std::atomic<std::size_t> next{0};

  auto worker = [&]() {
    std::vector<double> values(spec.axes.size(), 0.0);
    for (;;) {
      const std::size_t k = next.fetch_add(1);
      if (k >= total) return;
      decode(k, values);
      EngineConfig cfg = engine.config();
      cfg.sample_profile = false;
      for (std::size_t a = 0; a < spec.axes.size(); ++a) inputs[a][k] = values[a];
      // Everything that can throw has to stay inside the guard: an exception
      // escaping a worker thread would call std::terminate.
      try {
        for (std::size_t a = 0; a < spec.axes.size(); ++a)
          applyParameter(cfg, spec.axes[a].parameter, values[a]);
        const auto res = engine.runWith(cfg);
        for (std::size_t m = 0; m < spec.metrics.size(); ++m)
          outputs[m][k] = readMetric(res, spec.metrics[m]);
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

  SweepResult out;
  out.points = static_cast<int>(total);
  for (std::size_t a = 0; a < spec.axes.size(); ++a)
    out.table.addColumn(spec.axes[a].parameter, inputs[a], parameterUnits(spec.axes[a].parameter));
  for (std::size_t m = 0; m < spec.metrics.size(); ++m)
    out.table.addColumn(spec.metrics[m], outputs[m], metricUnits(spec.metrics[m]));
  out.table.addColumn("status", status);
  for (std::size_t k = 0; k < total; ++k) {
    if (status[k] == "failed") {
      ++out.failures;
      if (out.failure_messages.size() < 20) out.failure_messages.push_back(messages[k]);
    }
  }
  return out;
}

SweepSpec parseSweep(const ConfigNode& root) {
  SweepSpec spec;
  const auto s = root["sweep"];
  s.requireOnly({"axes", "metrics", "threads"});
  const auto axes = s["axes"];
  if (!axes.raw().IsSequence()) throw ConfigError("sweep: 'axes' must be a sequence");
  for (std::size_t i = 0; i < axes.raw().size(); ++i) {
    ConfigNode ax(axes.raw()[i], axes.path() + "[" + std::to_string(i) + "]", "sweep");
    ax.requireOnly({"parameter", "min", "max", "points", "spacing", "values"});
    const auto name = ax["parameter"].text();
    if (ax.has("values")) {
      SweepAxis a;
      a.parameter = name;
      a.values = ax["values"].numbers();
      spec.axes.push_back(a);
    } else {
      const auto spacing = ax.optional("spacing").text("linear");
      const double lo = ax["min"].number();
      const double hi = ax["max"].number();
      const int n = ax["points"].integer();
      if (spacing == "log" || spacing == "logarithmic") {
        spec.axes.push_back(SweepAxis::logarithmic(name, lo, hi, n));
      } else if (spacing == "linear") {
        spec.axes.push_back(SweepAxis::linear(name, lo, hi, n));
      } else {
        throw ConfigError("sweep axis '" + name + "': spacing must be 'linear' or 'log'");
      }
    }
  }
  spec.metrics = s["metrics"].strings();
  spec.threads = s.optional("threads").integer(0);
  return spec;
}

}  // namespace ignis
