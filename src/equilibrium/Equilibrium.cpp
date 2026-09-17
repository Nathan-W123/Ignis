// SPDX-License-Identifier: MIT
#include "ignis/equilibrium/Equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace ignis {

std::string toString(EquilibriumProblem p) {
  switch (p) {
    case EquilibriumProblem::kTP: return "TP";
    case EquilibriumProblem::kHP: return "HP";
    case EquilibriumProblem::kSP: return "SP";
    case EquilibriumProblem::kUV: return "UV";
    case EquilibriumProblem::kTV: return "TV";
  }
  return "?";
}

/// Active subset of elements and species for one problem instance.
///
/// Elements with b_i == 0 carry no atoms and would make the iteration matrix
/// singular; every species containing such an element is identically absent.
/// Both are removed up front and restored (as exact zeros) afterwards.
struct EquilibriumSolver::Reduced {
  std::vector<int> element_map;   // reduced index -> full element index
  std::vector<int> species_map;   // reduced index -> full species index
  Eigen::MatrixXd a;              // reduced element matrix (E' x N')
  Eigen::VectorXd b;              // reduced element moles per kg
  Eigen::VectorXd mw;             // reduced molar masses, kg/mol
  const SpeciesDatabase* db = nullptr;

  Eigen::Index numElements() const { return a.rows(); }
  Eigen::Index numSpecies() const { return a.cols(); }
  const Species& species(Eigen::Index j) const {
    return (*db)[static_cast<std::size_t>(species_map[static_cast<std::size_t>(j)])];
  }
  /// Expand a reduced mole-number vector to the full species list.
  Eigen::VectorXd expand(const Eigen::VectorXd& n_red) const {
    Eigen::VectorXd full = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db->size()));
    for (Eigen::Index j = 0; j < n_red.size(); ++j)
      full(species_map[static_cast<std::size_t>(j)]) = n_red(j);
    return full;
  }
};

EquilibriumSolver::EquilibriumSolver(const SpeciesDatabase& db, EquilibriumOptions opts)
    : db_(&db), mix_(db), opts_(std::move(opts)) {
  if (opts_.max_iterations < 1) throw ConfigError("equilibrium: max_iterations must be >= 1");
  if (!(opts_.species_tolerance > 0.0)) throw ConfigError("equilibrium: species_tolerance must be > 0");
}

double EquilibriumSolver::elementMass(const Eigen::VectorXd& b) const {
  double m = 0.0;
  for (Eigen::Index i = 0; i < b.size(); ++i)
    m += b(i) * db_->atomicWeight(db_->elements()[static_cast<std::size_t>(i)]);
  return m;
}

/// Build the reduced problem, dropping empty elements and the species that
/// depend on them.
void EquilibriumSolver::buildReduced(const SpeciesDatabase& db, const Eigen::VectorXd& b_full,
                  EquilibriumSolver::Reduced& red) {
  const Eigen::Index E = static_cast<Eigen::Index>(db.numElements());
  const Eigen::Index N = static_cast<Eigen::Index>(db.size());
  if (b_full.size() != E)
    throw ConfigError("equilibrium: element vector has length " +
                      std::to_string(b_full.size()) + ", expected " + std::to_string(E));
  if ((b_full.array() < 0.0).any())
    throw ConfigError("equilibrium: negative element mole number supplied");

  const double b_max = b_full.maxCoeff();
  if (!(b_max > 0.0)) throw ConfigError("equilibrium: element vector is all zeros");
  const double b_cut = 1.0e-14 * b_max;

  red.db = &db;
  red.element_map.clear();
  red.species_map.clear();
  std::vector<bool> element_active(static_cast<std::size_t>(E), false);
  for (Eigen::Index i = 0; i < E; ++i) {
    if (b_full(i) > b_cut) {
      element_active[static_cast<std::size_t>(i)] = true;
      red.element_map.push_back(static_cast<int>(i));
    }
  }
  const auto& a_full = db.elementMatrix();
  for (Eigen::Index j = 0; j < N; ++j) {
    bool ok = true;
    for (Eigen::Index i = 0; i < E; ++i) {
      if (a_full(i, j) > 0 && !element_active[static_cast<std::size_t>(i)]) { ok = false; break; }
    }
    if (ok) red.species_map.push_back(static_cast<int>(j));
  }
  if (red.species_map.empty())
    throw InfeasibleError("equilibrium: no species can be formed from the supplied elements");

  const Eigen::Index Er = static_cast<Eigen::Index>(red.element_map.size());
  const Eigen::Index Nr = static_cast<Eigen::Index>(red.species_map.size());
  red.a.setZero(Er, Nr);
  red.b.setZero(Er);
  red.mw.setZero(Nr);
  for (Eigen::Index i = 0; i < Er; ++i) {
    red.b(i) = b_full(red.element_map[static_cast<std::size_t>(i)]);
    for (Eigen::Index j = 0; j < Nr; ++j)
      red.a(i, j) = a_full(red.element_map[static_cast<std::size_t>(i)],
                           red.species_map[static_cast<std::size_t>(j)]);
  }
  for (Eigen::Index j = 0; j < Nr; ++j)
    red.mw(j) = db[static_cast<std::size_t>(red.species_map[static_cast<std::size_t>(j)])].molarMass();

  // Every active element must appear in at least one retained species,
  // otherwise its constraint row is identically zero and unsatisfiable.
  for (Eigen::Index i = 0; i < Er; ++i) {
    if (red.a.row(i).maxCoeff() <= 0.0) {
      throw InfeasibleError(
          "equilibrium: element '" +
          db.elements()[static_cast<std::size_t>(red.element_map[static_cast<std::size_t>(i)])] +
          "' is present in the reactants but in no available species");
    }
  }
}

EquilibriumResult EquilibriumSolver::attempt(const Reduced& red, EquilibriumProblem problem,
                                             double p, double target, double T_start,
                                             const Eigen::VectorXd& n_start) const {
  const Eigen::Index E = red.numElements();
  const Eigen::Index N = red.numSpecies();
  const bool solve_T = (problem != EquilibriumProblem::kTP);
  const Eigen::Index dim = E + 1 + (solve_T ? 1 : 0);
  const double R = constants::R_universal;
  const double T_lo = red.db->tMinCommon();
  const double T_hi = red.db->tMaxCommon();

  Eigen::VectorXd n = n_start;
  double T = T_start;
  double n_tot = n.sum();

  Eigen::VectorXd h_RT(N), s_R(N), cp_R(N), mu_RT(N);
  Eigen::MatrixXd A(dim, dim);
  Eigen::VectorXd rhs(dim), sol(dim), dln_n(N);

  EquilibriumDiagnostics diag;
  const double lnp = std::log(p / constants::p_reference);

  for (int iter = 0; iter < opts_.max_iterations; ++iter) {
    if (T < T_lo || T > T_hi) {
      std::ostringstream os;
      os << "temperature iterate " << T << " K left the polynomial range [" << T_lo << ", "
         << T_hi << "] K";
      throw ConvergenceError(os.str());
    }
    n_tot = n.sum();
    if (!(n_tot > 0.0)) throw ConvergenceError("total moles collapsed to zero");

    // --- species properties at the current temperature -----------------
    for (Eigen::Index j = 0; j < N; ++j) {
      double c, hh, ss;
      red.species(j).reduced(T, c, hh, ss);
      cp_R(j) = c;
      h_RT(j) = hh;
      // Entropy including the ideal-gas mixing and pressure terms.
      s_R(j) = ss - std::log(std::max(n(j), opts_.trace_floor) / n_tot) - lnp;
      mu_RT(j) = h_RT(j) - s_R(j);
    }

    // Weights: species far below the floor cannot influence the matrix.
    Eigen::VectorXd w = n;
    for (Eigen::Index j = 0; j < N; ++j) if (w(j) < opts_.trace_floor) w(j) = 0.0;

    // --- assemble the iteration matrix (RP-1311 Eqs. 2.24-2.26) --------
    A.setZero();
    rhs.setZero();
    const Eigen::VectorXd bw = red.a * w;                 // sum_j a_ij n_j
    A.topLeftCorner(E, E) = red.a * w.asDiagonal() * red.a.transpose();
    A.block(0, E, E, 1) = bw;
    A.block(E, 0, 1, E) = bw.transpose();
    A(E, E) = w.sum() - n_tot;
    rhs.head(E) = red.b - bw + red.a * (w.cwiseProduct(mu_RT));
    rhs(E) = n_tot - w.sum() + w.dot(mu_RT);

    if (solve_T) {
      Eigen::VectorXd key(N);
      if (problem == EquilibriumProblem::kHP) key = h_RT;
      else                                    key = s_R;      // kSP
      const Eigen::VectorXd wk = w.cwiseProduct(key);
      // Temperature column: derivative of Delta ln n_j wrt Delta ln T is h_RT.
      const Eigen::VectorXd wh = w.cwiseProduct(h_RT);
      A.block(0, E + 1, E, 1) = red.a * wh;
      A(E, E + 1) = wh.sum();
      // Energy / entropy row.
      A.block(E + 1, 0, 1, E) = (red.a * wk).transpose();
      A(E + 1, E) = wk.sum();
      A(E + 1, E + 1) = w.dot(cp_R) + wk.dot(h_RT);

      double closure;
      if (problem == EquilibriumProblem::kHP) {
        const double h = R * T * n.dot(h_RT);                // J/kg
        closure = (target - h) / (R * T);
      } else {
        const double s = R * n.dot(s_R);                     // J/(kg K)
        closure = (target - s) / R;
      }
      rhs(E + 1) = closure + wk.dot(mu_RT);
    }

    // --- solve ---------------------------------------------------------
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    if (!lu.isInvertible()) {
      std::ostringstream os;
      os << "iteration matrix is singular at iteration " << iter << " (T = " << T
         << " K, p = " << p << " Pa)";
      throw ConvergenceError(os.str());
    }
    sol = lu.solve(rhs);
    const Eigen::VectorXd pi = sol.head(E);
    const double dln_ntot = sol(E);
    const double dln_T = solve_T ? sol(E + 1) : 0.0;
    if (!sol.allFinite()) throw ConvergenceError("non-finite Newton correction");

    for (Eigen::Index j = 0; j < N; ++j)
      dln_n(j) = red.a.col(j).dot(pi) + dln_ntot + h_RT(j) * dln_T - mu_RT(j);

    // --- step-size control (RP-1311 Eqs. 3.1-3.3) ----------------------
    double worst = std::max(5.0 * std::abs(dln_ntot), 5.0 * std::abs(dln_T));
    for (Eigen::Index j = 0; j < N; ++j)
      if (n(j) / n_tot > 1.0e-8) worst = std::max(worst, std::abs(dln_n(j)));
    double lambda = (worst > 0.0) ? std::min(1.0, 2.0 / worst) : 1.0;
    for (Eigen::Index j = 0; j < N; ++j) {
      const double frac = n(j) / n_tot;
      if (frac <= 1.0e-8 && dln_n(j) > 0.0) {
        const double denom = std::abs(dln_n(j) - dln_ntot);
        if (denom > 0.0) {
          const double l2 = std::abs(std::log(1.0e-4 * n_tot / std::max(n(j), opts_.trace_floor))) / denom;
          lambda = std::min(lambda, l2);
        }
      }
    }
    lambda = std::max(lambda, 1.0e-6);

    // --- convergence test (on the *unlimited* Newton correction) -------
    double species_corr = 0.0;
    for (Eigen::Index j = 0; j < N; ++j)
      species_corr = std::max(species_corr, n(j) * std::abs(dln_n(j)) / n_tot);
    if (opts_.record_history) diag.history.push_back(species_corr);

    const bool done = species_corr < opts_.species_tolerance &&
                      std::abs(dln_ntot) < opts_.moles_tolerance &&
                      std::abs(dln_T) < opts_.temperature_tolerance;

    // --- apply the step ------------------------------------------------
    for (Eigen::Index j = 0; j < N; ++j) {
      const double ln_new = std::log(std::max(n(j), std::numeric_limits<double>::min())) +
                            lambda * dln_n(j);
      n(j) = (ln_new < -700.0) ? 0.0 : std::exp(ln_new);
      if (n(j) < opts_.trace_floor) n(j) = opts_.trace_floor;
    }
    if (solve_T) T *= std::exp(lambda * dln_T);

    if (done) {
      diag.converged = true;
      diag.iterations = iter + 1;
      diag.final_species_correction = species_corr;
      diag.final_moles_correction = std::abs(dln_ntot);
      diag.final_temperature_correction = std::abs(dln_T);
      break;
    }
    diag.iterations = iter + 1;
    diag.final_species_correction = species_corr;
    diag.final_moles_correction = std::abs(dln_ntot);
    diag.final_temperature_correction = std::abs(dln_T);
  }

  if (!diag.converged) {
    std::ostringstream os;
    os << toString(problem) << " equilibrium did not converge in " << opts_.max_iterations
       << " iterations (species correction " << std::scientific << std::setprecision(3)
       << diag.final_species_correction << ", |dln n| " << diag.final_moles_correction
       << ", |dln T| " << diag.final_temperature_correction << ", T = " << std::fixed
       << std::setprecision(2) << T << " K, p = " << p << " Pa)";
    throw ConvergenceError(os.str());
  }

  // --- residual audit ---------------------------------------------------
  n_tot = n.sum();
  const Eigen::VectorXd b_achieved = red.a * n;
  diag.element_residual_abs = (b_achieved - red.b).cwiseAbs().maxCoeff();
  diag.element_residual_rel = diag.element_residual_abs / red.b.maxCoeff();
  diag.mass_residual = n.dot(red.mw) - 1.0;

  // Recompute mu and pi at the converged point for the optimality residual.
  for (Eigen::Index j = 0; j < N; ++j) {
    double c, hh, ss;
    red.species(j).reduced(T, c, hh, ss);
    cp_R(j) = c;
    h_RT(j) = hh;
    s_R(j) = ss - std::log(std::max(n(j), opts_.trace_floor) / n_tot) - lnp;
    mu_RT(j) = h_RT(j) - s_R(j);
  }
  Eigen::VectorXd w = n;
  for (Eigen::Index j = 0; j < N; ++j) if (w(j) < opts_.trace_floor) w(j) = 0.0;
  // Least-squares recovery of pi from the stationarity conditions, weighted by
  // mole number so that trace species do not dominate.
  Eigen::MatrixXd M = red.a * w.asDiagonal() * red.a.transpose();
  Eigen::VectorXd rhs_pi = red.a * (w.cwiseProduct(mu_RT));
  Eigen::VectorXd pi = M.ldlt().solve(rhs_pi);

  diag.gibbs_residual = 0.0;
  for (Eigen::Index j = 0; j < N; ++j) {
    if (n(j) / n_tot < opts_.optimality_check_fraction) continue;
    diag.gibbs_residual =
        std::max(diag.gibbs_residual, std::abs(mu_RT(j) - red.a.col(j).dot(pi)));
  }

  EquilibriumResult res;
  res.pi = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db_->numElements()));
  for (Eigen::Index i = 0; i < E; ++i)
    res.pi(red.element_map[static_cast<std::size_t>(i)]) = pi(i);

  const Eigen::VectorXd n_full = red.expand(n);
  res.state = mix_.frozenState(n_full, T, p);
  res.state.model = CompositionModel::kEquilibrium;
  computeDerivatives(red, n, T, p, res.state);

  switch (problem) {
    case EquilibriumProblem::kHP:
      diag.state_residual = std::abs(res.state.h - target) / std::max(1.0, std::abs(target));
      break;
    case EquilibriumProblem::kSP:
      diag.state_residual = std::abs(res.state.s - target) / std::max(1.0, std::abs(target));
      break;
    default:
      diag.state_residual = 0.0;
      break;
  }

  if (diag.element_residual_rel > opts_.element_tolerance) {
    std::ostringstream os;
    os << toString(problem) << " equilibrium converged but element balance is off by "
       << std::scientific << diag.element_residual_rel << " (tolerance "
       << opts_.element_tolerance << ")";
    throw ConvergenceError(os.str());
  }
  if (diag.gibbs_residual > opts_.optimality_tolerance) {
    std::ostringstream os;
    os << toString(problem) << " equilibrium converged but Gibbs optimality residual is "
       << std::scientific << diag.gibbs_residual << " (tolerance "
       << opts_.optimality_tolerance << ")";
    throw ConvergenceError(os.str());
  }

  std::ostringstream msg;
  msg << toString(problem) << " converged in " << diag.iterations << " iterations";
  diag.message = msg.str();
  res.diagnostics = std::move(diag);
  return res;
}

void EquilibriumSolver::computeDerivatives(const Reduced& red, const Eigen::VectorXd& n_red,
                                           double T, double p, GasState& st) const {
  // Auxiliary systems from RP-1311 section 2.4.  At the converged point
  //   ln n_j = ln n - ln(p/p0) - g0_j/(R T) + sum_i a_ij pi_i,
  // so differentiating with respect to ln T (at fixed p) and ln p (at fixed T)
  // and imposing the element and total-mole constraints gives two linear
  // systems that share the matrix of the last Newton step.
  const Eigen::Index E = red.numElements();
  const Eigen::Index N = red.numSpecies();
  const double R = constants::R_universal;

  Eigen::VectorXd h_RT(N), cp_R(N), w = n_red;
  for (Eigen::Index j = 0; j < N; ++j) {
    double c, hh, ss;
    red.species(j).reduced(T, c, hh, ss);
    cp_R(j) = c;
    h_RT(j) = hh;
    if (w(j) < opts_.trace_floor) w(j) = 0.0;
  }
  const double n_tot = n_red.sum();

  Eigen::MatrixXd A(E + 1, E + 1);
  A.setZero();
  const Eigen::VectorXd bw = red.a * w;
  A.topLeftCorner(E, E) = red.a * w.asDiagonal() * red.a.transpose();
  A.block(0, E, E, 1) = bw;
  A.block(E, 0, 1, E) = bw.transpose();
  A(E, E) = w.sum() - n_tot;

  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
  if (!lu.isInvertible()) {
    // Fall back to the frozen derivative values rather than inventing numbers.
    st.cp_eff = st.cp_frozen;
    st.cv_eff = st.cv_frozen;
    st.gamma_eff = st.gamma_frozen;
    st.dlnV_dlnT_p = 1.0;
    st.dlnV_dlnp_T = -1.0;
    st.gamma_s = st.gamma_frozen;
    st.a = std::sqrt(st.gamma_s * p / st.rho);
    return;
  }

  // (d/d ln T)_p
  Eigen::VectorXd rhs_T(E + 1);
  const Eigen::VectorXd wh = w.cwiseProduct(h_RT);
  rhs_T.head(E) = -(red.a * wh);
  rhs_T(E) = -wh.sum();
  const Eigen::VectorXd sol_T = lu.solve(rhs_T);
  const double dlnn_dlnT = sol_T(E);
  Eigen::VectorXd dlnnj_dlnT(N);
  for (Eigen::Index j = 0; j < N; ++j)
    dlnnj_dlnT(j) = red.a.col(j).dot(sol_T.head(E)) + dlnn_dlnT + h_RT(j);

  // (d/d ln p)_T
  Eigen::VectorXd rhs_p(E + 1);
  rhs_p.head(E) = bw;
  rhs_p(E) = n_tot;
  const Eigen::VectorXd sol_p = lu.solve(rhs_p);
  const double dlnn_dlnp = sol_p(E);

  // Specific heat with composition shift:  cp = sum n_j cp_j + (1/T) sum n_j h_j dln n_j/dln T
  double cp = 0.0, shift = 0.0;
  for (Eigen::Index j = 0; j < N; ++j) {
    cp += n_red(j) * cp_R(j) * R;
    shift += n_red(j) * h_RT(j) * R * dlnnj_dlnT(j);
  }
  st.cp_eff = cp + shift;
  st.dlnV_dlnT_p = 1.0 + dlnn_dlnT;
  st.dlnV_dlnp_T = -1.0 + dlnn_dlnp;

  if (!(st.dlnV_dlnp_T < 0.0) || !(st.cp_eff > 0.0)) {
    st.cp_eff = st.cp_frozen;
    st.cv_eff = st.cv_frozen;
    st.gamma_eff = st.gamma_frozen;
    st.dlnV_dlnT_p = 1.0;
    st.dlnV_dlnp_T = -1.0;
    st.gamma_s = st.gamma_frozen;
    st.a = std::sqrt(st.gamma_s * p / st.rho);
    return;
  }

  // cp - cv = -(p v / T) (dlnV/dlnT)^2 / (dlnV/dlnp)
  const double v = st.v;
  st.cv_eff = st.cp_eff + (p * v / T) * st.dlnV_dlnT_p * st.dlnV_dlnT_p / st.dlnV_dlnp_T;
  st.gamma_eff = st.cp_eff / st.cv_eff;
  st.gamma_s = -st.gamma_eff / st.dlnV_dlnp_T;
  st.a = std::sqrt(st.gamma_s * p * v);
}

EquilibriumResult EquilibriumSolver::solveConstantPressure(const Eigen::VectorXd& b,
                                                           EquilibriumProblem problem, double p,
                                                           double target, double T_fixed_or_guess,
                                                           const Eigen::VectorXd* n_guess) const {
  if (!(p > 0.0)) throw ConfigError("equilibrium: pressure must be positive, got " + std::to_string(p));
  Reduced red;
  buildReduced(*db_, b, red);

  const Eigen::Index N = red.numSpecies();
  const double b_sum = red.b.sum();

  std::vector<double> temperatures;
  if (problem == EquilibriumProblem::kTP) {
    temperatures.push_back(T_fixed_or_guess);
  } else {
    temperatures.push_back(T_fixed_or_guess > 0.0 ? T_fixed_or_guess : opts_.initial_temperature);
    for (double t : opts_.restart_temperatures) temperatures.push_back(t);
  }

  std::vector<Eigen::VectorXd> starts;
  if (n_guess != nullptr) {
    Eigen::VectorXd ng(N);
    for (Eigen::Index j = 0; j < N; ++j)
      ng(j) = std::max((*n_guess)(red.species_map[static_cast<std::size_t>(j)]),
                       opts_.trace_floor);
    starts.push_back(ng);
  }
  starts.push_back(Eigen::VectorXd::Constant(N, 0.5 * b_sum / static_cast<double>(N)));
  starts.push_back(Eigen::VectorXd::Constant(N, 0.1 * b_sum / static_cast<double>(N)));
  starts.push_back(Eigen::VectorXd::Constant(N, 2.0 * b_sum / static_cast<double>(N)));

  std::ostringstream failures;
  int restarts = 0;
  for (std::size_t ti = 0; ti < temperatures.size(); ++ti) {
    for (std::size_t si = 0; si < starts.size(); ++si) {
      const double T0 = std::min(std::max(temperatures[ti], db_->tMinCommon() + 1.0),
                                 db_->tMaxCommon() - 1.0);
      try {
        EquilibriumResult r = attempt(red, problem, p, target, T0, starts[si]);
        r.diagnostics.restarts = restarts;
        return r;
      } catch (const IgnisError& e) {
        ++restarts;
        if (restarts <= 3) failures << "\n  attempt " << restarts << " (T0=" << T0 << " K): " << e.what();
      }
    }
  }
  std::ostringstream os;
  os << toString(problem) << " equilibrium failed after " << restarts << " attempts at p = " << p
     << " Pa. Element vector (mol/kg):";
  for (Eigen::Index i = 0; i < red.numElements(); ++i)
    os << " " << db_->elements()[static_cast<std::size_t>(red.element_map[static_cast<std::size_t>(i)])]
       << "=" << red.b(i);
  os << failures.str();
  throw ConvergenceError(os.str());
}

EquilibriumResult EquilibriumSolver::tp(const Eigen::VectorXd& b, double T, double p,
                                        const Eigen::VectorXd* n_guess) const {
  if (!(T > 0.0)) throw ConfigError("equilibrium: temperature must be positive");
  mix_.requireInRange(T);
  return solveConstantPressure(b, EquilibriumProblem::kTP, p, 0.0, T, n_guess);
}

EquilibriumResult EquilibriumSolver::hp(const Eigen::VectorXd& b, double h_target, double p,
                                        double T_guess, const Eigen::VectorXd* n_guess) const {
  return solveConstantPressure(b, EquilibriumProblem::kHP, p, h_target, T_guess, n_guess);
}

EquilibriumResult EquilibriumSolver::sp(const Eigen::VectorXd& b, double s_target, double p,
                                        double T_guess, const Eigen::VectorXd* n_guess) const {
  return solveConstantPressure(b, EquilibriumProblem::kSP, p, s_target, T_guess, n_guess);
}

EquilibriumResult EquilibriumSolver::uv(const Eigen::VectorXd& b, double u_target, double v_target,
                                        double T_guess) const {
  // Outer secant iteration on T.  For each trial T the pressure follows from
  // p = n R T / v where n is itself a function of the equilibrium solve, so an
  // inner fixed-point on p is run first.
  if (!(v_target > 0.0)) throw ConfigError("equilibrium (UV): specific volume must be positive");
  double T = (T_guess > 0.0) ? T_guess : opts_.initial_temperature;
  T = std::min(std::max(T, db_->tMinCommon() + 1.0), db_->tMaxCommon() - 1.0);

  EquilibriumResult r;
  double T_prev = T, f_prev = 0.0;
  int outer = 0;
  for (; outer < opts_.max_outer_iterations; ++outer) {
    // inner: converge pressure at this temperature
    double p = constants::R_universal * T / (v_target * 0.025);  // crude first guess, M ~ 25 g/mol
    for (int k = 0; k < 60; ++k) {
      r = tp(b, T, p);
      const double p_new = r.state.n_total * constants::R_universal * T / v_target;
      const double rel = std::abs(p_new - p) / p_new;
      p = p_new;
      if (rel < 1.0e-13) break;
      if (k == 59)
        throw ConvergenceError("equilibrium (UV): inner pressure fixed point did not converge");
    }
    const double f = r.state.u - u_target;
    const double scale = std::max(1.0, std::abs(u_target));
    if (std::abs(f) < 1.0e-11 * scale) {
      r.diagnostics.outer_iterations = outer + 1;
      r.diagnostics.state_residual = std::abs(f) / scale;
      return r;
    }
    double dT;
    if (outer == 0) {
      const double cv = std::max(r.state.cv_eff, 1.0);
      dT = -f / cv;
    } else {
      const double denom = (f - f_prev);
      dT = (std::abs(denom) > 0.0) ? -f * (T - T_prev) / denom : -f / std::max(r.state.cv_eff, 1.0);
    }
    dT = std::max(-500.0, std::min(500.0, dT));
    T_prev = T;
    f_prev = f;
    T = std::min(std::max(T + dT, db_->tMinCommon() + 1.0), db_->tMaxCommon() - 1.0);
  }
  throw ConvergenceError("equilibrium (UV): outer temperature iteration did not converge in " +
                         std::to_string(opts_.max_outer_iterations) + " steps");
}

EquilibriumResult EquilibriumSolver::tv(const Eigen::VectorXd& b, double T, double v_target,
                                        double p_guess) const {
  if (!(v_target > 0.0)) throw ConfigError("equilibrium (TV): specific volume must be positive");
  double p = (p_guess > 0.0) ? p_guess : constants::R_universal * T / (v_target * 0.025);
  EquilibriumResult r;
  for (int k = 0; k < opts_.max_outer_iterations; ++k) {
    r = tp(b, T, p);
    const double p_new = r.state.n_total * constants::R_universal * T / v_target;
    const double rel = std::abs(p_new - p) / p_new;
    p = p_new;
    if (rel < 1.0e-13) {
      r.diagnostics.outer_iterations = k + 1;
      return r;
    }
  }
  throw ConvergenceError("equilibrium (TV): pressure fixed point did not converge");
}

}  // namespace ignis
