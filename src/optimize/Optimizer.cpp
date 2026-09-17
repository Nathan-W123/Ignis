// SPDX-License-Identifier: MIT
#include "ignis/optimize/Optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>

#include "ignis/optimize/Parameters.hpp"

namespace ignis {

double DesignConstraint::violation(double value) const {
  const double s = (scale > 0.0) ? scale : std::max(1.0, std::abs(bound));
  if (op == "<=") return (value - bound) / s;
  if (op == ">=") return (bound - value) / s;
  throw ConfigError("constraint on '" + metric + "': operator must be '<=' or '>='");
}

namespace {

/// Evaluation of one design point.
struct Evaluation {
  bool ok = false;
  double objective = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> constraints;
  std::vector<double> violations;
  double max_violation = 0.0;
  SteadyEngineResult result;
  std::string message;
};

Evaluation evaluate(const SteadyEngine& engine, const OptimizationSpec& spec,
                    const std::vector<double>& x) {
  Evaluation e;
  e.constraints.assign(spec.constraints.size(), std::numeric_limits<double>::quiet_NaN());
  e.violations.assign(spec.constraints.size(), 0.0);
  EngineConfig cfg = engine.config();
  cfg.sample_profile = false;
  for (std::size_t i = 0; i < spec.variables.size(); ++i)
    applyParameter(cfg, spec.variables[i].parameter, x[i]);
  try {
    e.result = engine.runWith(cfg);
    e.objective = readMetric(e.result, spec.objective);
    e.max_violation = 0.0;
    for (std::size_t c = 0; c < spec.constraints.size(); ++c) {
      e.constraints[c] = readMetric(e.result, spec.constraints[c].metric);
      e.violations[c] = spec.constraints[c].violation(e.constraints[c]);
      e.max_violation = std::max(e.max_violation, e.violations[c]);
    }
    e.ok = std::isfinite(e.objective);
    if (!e.ok) e.message = "objective metric was not finite";
  } catch (const std::exception& ex) {
    e.ok = false;
    e.message = ex.what();
  }
  return e;
}

/// Augmented-Lagrangian merit value (always minimised).
double merit(const OptimizationSpec& spec, const Evaluation& e, const std::vector<double>& lambda,
             double mu, double fail_penalty) {
  if (!e.ok) return fail_penalty;
  double m = (spec.maximize ? -e.objective : e.objective);
  for (std::size_t c = 0; c < spec.constraints.size(); ++c) {
    const double t = lambda[c] + mu * e.violations[c];
    const double pos = std::max(0.0, t);
    m += (pos * pos - lambda[c] * lambda[c]) / (2.0 * mu);
  }
  return m;
}

/// Latin hypercube sample in the unit cube.
std::vector<std::vector<double>> latinHypercube(int n, int dim, std::mt19937_64& rng) {
  std::vector<std::vector<double>> out(static_cast<std::size_t>(n),
                                       std::vector<double>(static_cast<std::size_t>(dim)));
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  for (int d = 0; d < dim; ++d) {
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int i = 0; i < n; ++i)
      out[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] =
          (perm[static_cast<std::size_t>(i)] + uni(rng)) / n;
  }
  return out;
}

}  // namespace

OptimizationResult optimize(const SteadyEngine& engine, const OptimizationSpec& spec) {
  if (spec.variables.empty()) throw ConfigError("optimization: no design variables");
  if (spec.objective.empty()) throw ConfigError("optimization: no objective metric");
  (void)metricUnits(spec.objective);
  for (const auto& v : spec.variables) {
    (void)parameterUnits(v.parameter);
    if (!(v.max > v.min))
      throw ConfigError("optimization: variable '" + v.parameter + "' has max <= min");
  }
  for (const auto& c : spec.constraints) {
    (void)metricUnits(c.metric);
    if (c.op != "<=" && c.op != ">=")
      throw ConfigError("optimization: constraint on '" + c.metric +
                        "' must use '<=' or '>='");
  }

  const std::size_t n = spec.variables.size();
  std::vector<double> lo(n), hi(n), x0(n);
  for (std::size_t i = 0; i < n; ++i) {
    lo[i] = spec.variables[i].min;
    hi[i] = spec.variables[i].max;
    const double s = spec.variables[i].start > 0.0
                         ? spec.variables[i].start
                         : readParameter(engine.config(), spec.variables[i].parameter);
    x0[i] = std::min(std::max(s, lo[i]), hi[i]);
  }
  auto toPhysical = [&](const std::vector<double>& u) {
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) {
      // Reflect at the bounds so the simplex can never leave the box.
      double t = u[i];
      for (int k = 0; k < 8 && (t < 0.0 || t > 1.0); ++k) {
        if (t < 0.0) t = -t;
        if (t > 1.0) t = 2.0 - t;
      }
      t = std::min(1.0, std::max(0.0, t));
      x[i] = lo[i] + t * (hi[i] - lo[i]);
    }
    return x;
  };
  auto toUnit = [&](const std::vector<double>& x) {
    std::vector<double> u(n);
    for (std::size_t i = 0; i < n; ++i) u[i] = (x[i] - lo[i]) / (hi[i] - lo[i]);
    return u;
  };

  OptimizationResult res;
  for (const auto& v : spec.variables) res.variable_names.push_back(v.parameter);

  // Dimension-adaptive Nelder-Mead coefficients (Gao & Han, 2012).
  const double dn = static_cast<double>(n);
  const double alpha = 1.0;
  const double beta = 1.0 + 2.0 / dn;
  const double gamma = 0.75 - 1.0 / (2.0 * dn);
  const double delta = 1.0 - 1.0 / dn;

  std::mt19937_64 rng(spec.seed);
  auto starts = latinHypercube(std::max(0, spec.starts - 1), static_cast<int>(n), rng);
  starts.insert(starts.begin(), toUnit(x0));

  double best_merit = std::numeric_limits<double>::infinity();
  Evaluation best_eval;
  std::vector<double> best_x = x0;
  bool have_best = false;
  int evaluations = 0;

  // A failed run must still produce a finite, large merit so the simplex can
  // walk out of the infeasible region.
  const double fail_penalty = 1.0e12;

  for (std::size_t s = 0; s < starts.size() && evaluations < spec.max_evaluations; ++s) {
    std::vector<double> lambda(spec.constraints.size(), 0.0);
    double mu = spec.initial_penalty;
    std::vector<double> u = starts[s];
    double prev_violation = std::numeric_limits<double>::infinity();

    for (int outer = 0; outer < spec.outer_iterations && evaluations < spec.max_evaluations;
         ++outer) {
      // --- build the initial simplex -----------------------------------
      std::vector<std::vector<double>> simplex;
      simplex.push_back(u);
      for (std::size_t i = 0; i < n; ++i) {
        auto v = u;
        v[i] += (v[i] < 0.5) ? 0.10 : -0.10;
        simplex.push_back(v);
      }
      std::vector<double> fvals(simplex.size());
      std::vector<Evaluation> evals(simplex.size());
      auto score = [&](const std::vector<double>& uu, Evaluation& e) {
        e = evaluate(engine, spec, toPhysical(uu));
        ++evaluations;
        if (!e.ok) ++res.failures;
        const double m = merit(spec, e, lambda, mu, fail_penalty);
        OptimizationStep st;
        st.evaluation = evaluations;
        st.outer = outer;
        st.objective = e.objective;
        st.merit = m;
        st.max_violation = e.max_violation;
        st.feasible = e.ok && e.max_violation <= 1e-9;
        st.x = toPhysical(uu);
        res.history.push_back(std::move(st));
        if (e.ok && e.max_violation <= 1e-9) {
          const double raw = spec.maximize ? -e.objective : e.objective;
          if (!have_best || raw < best_merit) {
            best_merit = raw;
            best_eval = e;
            best_x = toPhysical(uu);
            have_best = true;
          }
        }
        return m;
      };
      for (std::size_t i = 0; i < simplex.size(); ++i) fvals[i] = score(simplex[i], evals[i]);

      // --- Nelder-Mead --------------------------------------------------
      for (int it = 0; it < spec.max_evaluations && evaluations < spec.max_evaluations; ++it) {
        std::vector<std::size_t> order(simplex.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return fvals[a] < fvals[b]; });
        const std::size_t ib = order.front(), iw = order.back();
        const std::size_t is = order[order.size() - 2];

        double spread = 0.0;
        for (std::size_t i = 0; i < simplex.size(); ++i)
          for (std::size_t d = 0; d < n; ++d)
            spread = std::max(spread, std::abs(simplex[i][d] - simplex[ib][d]));
        if (spread < spec.simplex_tolerance) break;

        std::vector<double> centroid(n, 0.0);
        for (std::size_t i = 0; i < simplex.size(); ++i) {
          if (i == iw) continue;
          for (std::size_t d = 0; d < n; ++d) centroid[d] += simplex[i][d];
        }
        for (std::size_t d = 0; d < n; ++d) centroid[d] /= dn;

        auto combine = [&](double coeff) {
          std::vector<double> p(n);
          for (std::size_t d = 0; d < n; ++d)
            p[d] = centroid[d] + coeff * (centroid[d] - simplex[iw][d]);
          return p;
        };

        Evaluation er;
        auto xr = combine(alpha);
        const double fr = score(xr, er);
        if (fr < fvals[ib]) {
          Evaluation ee;
          auto xe = combine(beta);
          const double fe = score(xe, ee);
          if (fe < fr) { simplex[iw] = xe; fvals[iw] = fe; evals[iw] = ee; }
          else { simplex[iw] = xr; fvals[iw] = fr; evals[iw] = er; }
        } else if (fr < fvals[is]) {
          simplex[iw] = xr;
          fvals[iw] = fr;
          evals[iw] = er;
        } else {
          const bool outside = fr < fvals[iw];
          Evaluation ec;
          auto xc = combine(outside ? gamma : -gamma);
          const double fc = score(xc, ec);
          if (fc < (outside ? fr : fvals[iw])) {
            simplex[iw] = xc;
            fvals[iw] = fc;
            evals[iw] = ec;
          } else {
            for (std::size_t i = 0; i < simplex.size(); ++i) {
              if (i == ib) continue;
              for (std::size_t d = 0; d < n; ++d)
                simplex[i][d] = simplex[ib][d] + delta * (simplex[i][d] - simplex[ib][d]);
              fvals[i] = score(simplex[i], evals[i]);
            }
          }
        }
      }

      // --- augmented-Lagrangian update -----------------------------------
      std::size_t ib = 0;
      for (std::size_t i = 1; i < fvals.size(); ++i)
        if (fvals[i] < fvals[ib]) ib = i;
      u = simplex[ib];
      const Evaluation& e = evals[ib];
      if (spec.constraints.empty()) break;
      if (!e.ok) { mu *= spec.penalty_growth; continue; }
      for (std::size_t c = 0; c < spec.constraints.size(); ++c)
        lambda[c] = std::max(0.0, lambda[c] + mu * e.violations[c]);
      if (e.max_violation > 0.25 * prev_violation) mu *= spec.penalty_growth;
      prev_violation = std::max(0.0, e.max_violation);
      if (prev_violation < 1e-10) break;
    }
  }

  res.evaluations = evaluations;
  if (!have_best) {
    std::ostringstream os;
    os << "optimization: no feasible design was found in " << evaluations << " evaluations ("
       << res.failures << " analyses failed outright). Relax the constraints or widen the "
          "variable bounds.";
    throw InfeasibleError(os.str());
  }
  res.x = best_x;
  res.objective = best_eval.objective;
  res.constraint_values = best_eval.constraints;
  res.constraint_violations = best_eval.violations;
  res.feasible = true;
  res.best = best_eval.result;
  return res;
}

std::string OptimizationResult::summary(const OptimizationSpec& spec) const {
  std::ostringstream os;
  os << "optimization result\n"
     << "  objective            " << (spec.maximize ? "maximise " : "minimise ") << spec.objective
     << " = " << std::fixed << std::setprecision(6) << objective << " "
     << metricUnits(spec.objective) << "\n"
     << "  evaluations          " << evaluations << " (" << failures << " failed analyses)\n"
     << "  feasible             " << (feasible ? "yes" : "no") << "\n"
     << "  design variables\n";
  for (std::size_t i = 0; i < x.size(); ++i)
    os << "    " << std::left << std::setw(34) << variable_names[i] << std::right << std::setw(16)
       << x[i] << "  " << parameterUnits(variable_names[i]) << "\n";
  if (!spec.constraints.empty()) {
    os << "  constraints\n";
    for (std::size_t c = 0; c < spec.constraints.size(); ++c)
      os << "    " << std::left << std::setw(34) << spec.constraints[c].metric << std::right
         << std::setw(16) << constraint_values[c] << "  " << spec.constraints[c].op << " "
         << spec.constraints[c].bound << "  (normalised violation "
         << std::scientific << std::setprecision(2) << constraint_violations[c] << std::fixed
         << std::setprecision(6) << ")\n";
  }
  return os.str();
}

Table OptimizationResult::historyTable(const OptimizationSpec& spec) const {
  Table t("optimization_history");
  std::vector<double> ev, outer, obj, mer, viol, feas;
  std::vector<std::vector<double>> vars(variable_names.size());
  for (const auto& h : history) {
    ev.push_back(h.evaluation);
    outer.push_back(h.outer);
    obj.push_back(h.objective);
    mer.push_back(h.merit);
    viol.push_back(h.max_violation);
    feas.push_back(h.feasible ? 1.0 : 0.0);
    for (std::size_t i = 0; i < variable_names.size(); ++i) vars[i].push_back(h.x[i]);
  }
  t.addColumn("evaluation", ev, "-");
  t.addColumn("outer_iteration", outer, "-");
  for (std::size_t i = 0; i < variable_names.size(); ++i)
    t.addColumn(variable_names[i], vars[i], parameterUnits(variable_names[i]));
  t.addColumn(spec.objective, obj, metricUnits(spec.objective));
  t.addColumn("merit", mer, "-");
  t.addColumn("max_violation", viol, "-");
  t.addColumn("feasible", feas, "-");
  return t;
}

OptimizationSpec parseOptimization(const ConfigNode& root) {
  OptimizationSpec spec;
  const auto o = root["optimization"];
  o.requireOnly({"objective", "sense", "variables", "constraints", "starts", "max_evaluations",
                 "simplex_tolerance", "outer_iterations", "initial_penalty", "penalty_growth",
                 "seed", "threads"});
  spec.objective = o["objective"].text();
  const auto sense = o.optional("sense").text("maximize");
  if (sense != "maximize" && sense != "minimize")
    throw ConfigError("optimization: 'sense' must be 'maximize' or 'minimize'");
  spec.maximize = (sense == "maximize");

  const auto vars = o["variables"];
  if (!vars.raw().IsSequence()) throw ConfigError("optimization: 'variables' must be a sequence");
  for (std::size_t i = 0; i < vars.raw().size(); ++i) {
    ConfigNode v(vars.raw()[i], vars.path() + "[" + std::to_string(i) + "]", "optimization");
    v.requireOnly({"parameter", "min", "max", "start"});
    DesignVariable dv;
    dv.parameter = v["parameter"].text();
    dv.min = v["min"].number();
    dv.max = v["max"].number();
    dv.start = v.optional("start").number(0.0);
    spec.variables.push_back(dv);
  }
  if (o.has("constraints")) {
    const auto cs = o["constraints"];
    for (std::size_t i = 0; i < cs.raw().size(); ++i) {
      ConfigNode c(cs.raw()[i], cs.path() + "[" + std::to_string(i) + "]", "optimization");
      c.requireOnly({"metric", "op", "bound", "scale"});
      DesignConstraint dc;
      dc.metric = c["metric"].text();
      dc.op = c.optional("op").text("<=");
      dc.bound = c["bound"].number();
      dc.scale = c.optional("scale").number(0.0);
      spec.constraints.push_back(dc);
    }
  }
  spec.starts = o.optional("starts").integer(4);
  spec.max_evaluations = o.optional("max_evaluations").integer(4000);
  spec.simplex_tolerance = o.optional("simplex_tolerance").number(1e-12, 1e-1, 1.0e-7);
  spec.outer_iterations = o.optional("outer_iterations").integer(8);
  spec.initial_penalty = o.optional("initial_penalty").number(1e-3, 1e9, 10.0);
  spec.penalty_growth = o.optional("penalty_growth").number(1.1, 100.0, 5.0);
  spec.seed = static_cast<unsigned>(o.optional("seed").integer(20260917));
  spec.threads = o.optional("threads").integer(0);
  return spec;
}

}  // namespace ignis
