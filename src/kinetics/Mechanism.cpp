// SPDX-License-Identifier: MIT
#include "ignis/kinetics/Mechanism.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <set>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/core/Version.hpp"

namespace ignis {

double ArrheniusRate::operator()(double T) const {
  double k = A;
  if (b != 0.0) k *= std::pow(T, b);
  if (Ea != 0.0) k *= std::exp(-Ea / (constants::R_universal * T));
  return k;
}

namespace {

ArrheniusRate readRate(const YAML::Node& n, const std::string& where) {
  if (!n || !n.IsMap()) throw ConfigError("mechanism: " + where + " needs {A, b, Ea}");
  ArrheniusRate r;
  r.A = n["A"].as<double>();
  r.b = n["b"].as<double>(0.0);
  r.Ea = n["Ea"].as<double>(0.0);
  if (!(r.A >= 0.0)) throw ConfigError("mechanism: negative pre-exponential factor in " + where);
  return r;
}

}  // namespace

std::string Mechanism::findDefault() {
  namespace fs = std::filesystem;
  std::vector<std::string> tried;
  auto probe = [&](const fs::path& p) {
    tried.push_back(p.string());
    std::error_code ec;
    return fs::exists(p, ec) && fs::is_regular_file(p, ec);
  };
  const char* env = std::getenv("IGNIS_DATA_DIR");
  if (env != nullptr) {
    const fs::path p = fs::path(env) / "kinetics" / "gri30_nozzle.yaml";
    if (probe(p)) return p.string();
  }
  {
    const fs::path p = fs::path(kDefaultDataDir) / "kinetics" / "gri30_nozzle.yaml";
    if (probe(p)) return p.string();
  }
  for (const char* rel : {"data/kinetics/gri30_nozzle.yaml", "../data/kinetics/gri30_nozzle.yaml",
                          "../../data/kinetics/gri30_nozzle.yaml"}) {
    if (probe(rel)) return std::string(rel);
  }
  std::ostringstream os;
  os << "cannot locate the nozzle reaction mechanism. Tried:";
  for (const auto& t : tried) os << "\n  " << t;
  throw ConfigError(os.str());
}

Mechanism Mechanism::loadYaml(const std::string& path, const SpeciesDatabase& db) {
  YAML::Node doc;
  try {
    doc = YAML::LoadFile(path);
  } catch (const YAML::Exception& e) {
    throw ConfigError("cannot read mechanism '" + path + "': " + e.what());
  }
  Mechanism m;
  m.db_ = &db;
  m.path_ = path;
  const auto reactions = doc["reactions"];
  if (!reactions || !reactions.IsSequence())
    throw ConfigError("mechanism '" + path + "': no 'reactions' list");
  std::set<int> active;
  for (const auto& node : reactions) {
    Reaction r;
    r.equation = node["equation"].as<std::string>("?");
    r.source_id = node["gri"].as<int>(0);
    const std::string where = "'" + r.equation + "'";
    const std::string type = node["type"].as<std::string>("elementary");
    if (type == "elementary") r.type = ReactionType::kElementary;
    else if (type == "three-body") r.type = ReactionType::kThreeBody;
    else if (type == "falloff") r.type = ReactionType::kFalloff;
    else throw ConfigError("mechanism: unknown reaction type '" + type + "' in " + where);

    // A reaction naming a species this gas cannot contain (CO in a
    // hydrogen-oxygen engine) cannot occur in it, and is left out.
    std::string missing;
    auto side = [&](const char* key, std::vector<int>& idx, std::vector<double>& nu) {
      const auto s = node[key];
      if (!s || !s.IsMap() || s.size() == 0)
        throw ConfigError(std::string("mechanism: ") + where + " has no " + key);
      for (const auto& kv : s) {
        const auto name = kv.first.as<std::string>();
        const int j = db.index(name);
        if (j < 0) {
          missing = name;
          continue;
        }
        idx.push_back(j);
        nu.push_back(kv.second.as<double>());
      }
    };
    side("reactants", r.reactants, r.reactant_nu);
    side("products", r.products, r.product_nu);
    if (!missing.empty()) {
      m.dropped_.push_back(r.equation + " (" + missing + " is not in the species set)");
      continue;
    }
    for (double v : r.product_nu) r.delta_nu += v;
    for (double v : r.reactant_nu) r.delta_nu -= v;

    if (r.type == ReactionType::kFalloff) {
      r.rate = readRate(node["high"], where + " high");
      r.low = readRate(node["low"], where + " low");
      if (const auto t = node["troe"]) {
        if (!t.IsSequence() || t.size() < 3)
          throw ConfigError("mechanism: Troe parameters of " + where + " need [a, T3, T1(, T2)]");
        r.troe = true;
        r.troe_a = t[0].as<double>();
        r.troe_t3 = t[1].as<double>();
        r.troe_t1 = t[2].as<double>();
        if (t.size() > 3) {
          r.troe_t2 = t[3].as<double>();
          r.troe_has_t2 = true;
        }
      }
    } else {
      r.rate = readRate(node["rate"], where);
    }

    if (r.type != ReactionType::kElementary) {
      r.default_efficiency = node["default_efficiency"].as<double>(1.0);
      bool any = r.default_efficiency > 0.0;
      if (const auto e = node["efficiencies"]) {
        for (const auto& kv : e) {
          const int j = db.index(kv.first.as<std::string>());
          if (j < 0) continue;  // absent from this gas: zero concentration
          r.efficiencies.emplace_back(j, kv.second.as<double>());
          if (kv.second.as<double>() > 0.0) any = true;
        }
      }
      if (!any) {
        m.dropped_.push_back(r.equation + " (its collider is not in the species set)");
        continue;
      }
    }
    for (const int j : r.reactants) active.insert(j);
    for (const int j : r.products) active.insert(j);
    m.reactions_.push_back(std::move(r));
  }
  if (m.reactions_.empty())
    throw ConfigError("mechanism '" + path + "': no reaction can occur among the species set");
  m.active_.assign(active.begin(), active.end());
  return m;
}

namespace {

/// k_f and K_c of one reaction at T, with the third-body concentration M.
/// g0 holds the standard Gibbs energies at T in database order (only the
/// reaction's own species are read), and ln_p0_rt is ln(p0 / R T): both are
/// the same for every reaction, so the caller evaluates them once.
void forwardAndEquilibrium(const Reaction& r, const std::vector<double>& g0, double T,
                           double ln_p0_rt, double M, double& kf, double& kc) {
  if (r.type == ReactionType::kFalloff) {
    const double k_inf = r.rate(T);
    const double k_0 = r.low(T);
    const double pr = (k_inf > 0.0) ? k_0 * M / k_inf : 0.0;
    if (!(pr > 0.0)) {
      kf = 0.0;
    } else {
      double F = 1.0;
      if (r.troe) {
        double fcent = (1.0 - r.troe_a) * std::exp(-T / r.troe_t3) + r.troe_a * std::exp(-T / r.troe_t1);
        if (r.troe_has_t2) fcent += std::exp(-r.troe_t2 / T);
        fcent = std::max(fcent, 1e-300);
        const double lf = std::log10(fcent);
        const double c = -0.4 - 0.67 * lf;
        const double n = 0.75 - 1.27 * lf;
        const double lp = std::log10(pr) + c;
        const double f1 = lp / (n - 0.14 * lp);
        F = std::pow(10.0, lf / (1.0 + f1 * f1));
      }
      kf = k_inf * (pr / (1.0 + pr)) * F;
    }
  } else {
    kf = r.rate(T);
  }
  // Detailed balance on the species data: K_c = exp(-dG0/RT) (p0/RT)^dnu.
  double dg = 0.0;
  for (std::size_t i = 0; i < r.products.size(); ++i)
    dg += r.product_nu[i] * g0[static_cast<std::size_t>(r.products[i])];
  for (std::size_t i = 0; i < r.reactants.size(); ++i)
    dg -= r.reactant_nu[i] * g0[static_cast<std::size_t>(r.reactants[i])];
  const double RT = constants::R_universal * T;
  kc = std::exp(-dg / RT + r.delta_nu * ln_p0_rt);
}

double lnReferenceConcentration(double T) {
  return std::log(constants::p_reference / (constants::R_universal * T));
}

double thirdBody(const Reaction& r, const Eigen::VectorXd& C, double total) {
  double M = r.default_efficiency * total;
  for (const auto& [j, eps] : r.efficiencies)
    M += (eps - r.default_efficiency) * std::max(C(j), 0.0);
  return M;
}

}  // namespace

void Mechanism::rateConstants(std::size_t idx, double T, double M, double& kf, double& kr) const {
  const auto& r = reactions_.at(idx);
  std::vector<double> g0(db_->size(), 0.0);
  for (const int j : r.reactants) g0[static_cast<std::size_t>(j)] = (*db_)[static_cast<std::size_t>(j)].g0(T);
  for (const int j : r.products) g0[static_cast<std::size_t>(j)] = (*db_)[static_cast<std::size_t>(j)].g0(T);
  double kc = 0.0;
  forwardAndEquilibrium(r, g0, T, lnReferenceConcentration(T), M, kf, kc);
  kr = kf / kc;
}

void Mechanism::productionRates(double T, const Eigen::VectorXd& C, Eigen::VectorXd& wdot,
                                double rate_multiplier) const {
  const auto& db = *db_;
  wdot.setZero(static_cast<Eigen::Index>(db.size()));
  if (rate_multiplier == 0.0) return;
  double total = 0.0;
  for (Eigen::Index j = 0; j < C.size(); ++j) total += std::max(C(j), 0.0);
  // Each species' Gibbs energy is needed by every reaction it takes part in;
  // evaluate it once.
  std::vector<double> g0(db.size(), 0.0);
  for (const int j : active_) g0[static_cast<std::size_t>(j)] = db[static_cast<std::size_t>(j)].g0(T);
  const double ln_p0_rt = lnReferenceConcentration(T);
  for (const auto& r : reactions_) {
    const double M = (r.type == ReactionType::kElementary) ? 1.0 : thirdBody(r, C, total);
    double kf = 0.0, kc = 1.0;
    forwardAndEquilibrium(r, g0, T, ln_p0_rt, M, kf, kc);
    if (kf == 0.0) continue;
    double fwd = kf, rev = kf / kc;
    for (std::size_t i = 0; i < r.reactants.size(); ++i)
      fwd *= std::pow(std::max(C(r.reactants[i]), 0.0), r.reactant_nu[i]);
    for (std::size_t i = 0; i < r.products.size(); ++i)
      rev *= std::pow(std::max(C(r.products[i]), 0.0), r.product_nu[i]);
    double q = fwd - rev;
    if (r.type == ReactionType::kThreeBody) q *= M;
    q *= rate_multiplier;
    for (std::size_t i = 0; i < r.reactants.size(); ++i) wdot(r.reactants[i]) -= r.reactant_nu[i] * q;
    for (std::size_t i = 0; i < r.products.size(); ++i) wdot(r.products[i]) += r.product_nu[i] * q;
  }
}

}  // namespace ignis
