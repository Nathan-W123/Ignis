// SPDX-License-Identifier: MIT
/// \file ignis_equilibrium.cpp
/// \brief Solve and report a chemical-equilibrium state.

#include <algorithm>
#include <iomanip>
#include <cstdio>
#include <iostream>

#include "AppCommon.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/io/Table.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/thermo/Transport.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"problem", "NAME", "equilibrium problem: hp (adiabatic), tp, sp", "hp"},
        {"temperature", "K", "temperature for the tp problem", "from the chamber solution"},
        {"min-fraction", "X", "smallest mole fraction printed", "1e-8"},
        {"no-transport", "", "skip the transport-property evaluation", ""},
        {"mr-sweep", "MIN:MAX:N", "also write the composition over a mixture-ratio sweep", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_equilibrium",
                          "Solve chemical equilibrium for the propellant mixture defined by a\n"
                          "configuration file and report composition, thermodynamic properties\n"
                          "and solver residuals.",
                          options, cli))
      return 0;

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);

    const std::string db_path = cfg.species_database.empty() ? findDefaultSpeciesDatabase()
                                                             : cfg.species_database;
    const auto full = SpeciesDatabase::loadYaml(db_path);
    const SpeciesDatabase db = cfg.species_include.empty()
                                   ? full.restrictToElements(cfg.elements)
                                   : full.subset(cfg.species_include);
    const EquilibriumSolver solver(db);
    const auto lib = PropellantLibrary::loadYaml(
        cfg.propellant_library.empty() ? findDefaultPropellantLibrary() : cfg.propellant_library);
    const auto& ox = lib.at(cfg.oxidizer);
    const auto& fu = lib.at(cfg.fuel);
    const PropellantMixture mix(
        ox, fu, cfg.mixture_ratio,
        cfg.oxidizer_temperature > 0.0 ? cfg.oxidizer_temperature : ox.reference_temperature,
        cfg.fuel_temperature > 0.0 ? cfg.fuel_temperature : fu.reference_temperature);

    const Eigen::VectorXd b = mix.elementMoles(db);
    const double h0 = mix.enthalpy(db);
    const std::string problem = flagString(cli, "problem", "hp");

    EquilibriumResult eq;
    if (problem == "hp") {
      eq = solver.hp(b, h0, cfg.chamber_pressure);
    } else if (problem == "tp") {
      const double T = flagDouble(cli, "temperature", 0.0);
      if (!(T > 0.0))
        throw ConfigError("--problem tp needs --temperature to be given in kelvin");
      eq = solver.tp(b, T, cfg.chamber_pressure);
    } else if (problem == "sp") {
      const auto ref = solver.hp(b, h0, cfg.chamber_pressure);
      eq = solver.sp(b, ref.state.s, cfg.chamber_pressure, ref.state.T);
    } else {
      throw ConfigError("--problem must be 'hp', 'tp' or 'sp', got '" + problem + "'");
    }

    const double min_fraction = flagDouble(cli, "min-fraction", 1.0e-8);
    const auto X = eq.state.moleFractions();
    const auto Y = eq.state.massFractions(db.molarMasses());

    std::vector<std::size_t> order(db.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) {
      return X(static_cast<Eigen::Index>(a)) > X(static_cast<Eigen::Index>(c));
    });

    if (!cli.quiet) {
      std::cout << versionBanner() << "\n\n"
                << "case            " << cfg.name << "\n"
                << "species data    " << db.provenance().path << "\n";
      for (const auto& s : db.provenance().sources) std::cout << "                " << s << "\n";
      std::cout << "problem         " << problem << "\n"
                << "propellants     " << cfg.oxidizer << " / " << cfg.fuel << "  O/F = "
                << std::fixed << std::setprecision(4) << cfg.mixture_ratio
                << "  (phi = " << mix.equivalenceRatio() << ")\n"
                << "reactant h      " << std::setprecision(1) << h0 * 1e-3 << " kJ/kg\n"
                << "pressure        " << std::setprecision(5) << cfg.chamber_pressure * 1e-6
                << " MPa\n"
                << "species in set  " << db.size() << " over elements";
      for (const auto& e : db.elements()) std::cout << " " << e;
      std::cout << "\n\n"
                << "temperature     " << std::setprecision(3) << eq.state.T << " K\n"
                << "molar mass      " << std::setprecision(5) << eq.state.M * 1e3 << " g/mol\n"
                << "density         " << eq.state.rho << " kg/m^3\n"
                << "enthalpy        " << std::setprecision(1) << eq.state.h * 1e-3 << " kJ/kg\n"
                << "entropy         " << std::setprecision(3) << eq.state.s << " J/(kg K)\n"
                << "cp frozen/eq    " << std::setprecision(2) << eq.state.cp_frozen << " / "
                << eq.state.cp_eff << " J/(kg K)\n"
                << "cv frozen/eq    " << eq.state.cv_frozen << " / " << eq.state.cv_eff
                << " J/(kg K)\n"
                << "gamma frozen    " << std::setprecision(6) << eq.state.gamma_frozen << "\n"
                << "gamma_s         " << eq.state.gamma_s << "\n"
                << "(dlnV/dlnT)_p   " << eq.state.dlnV_dlnT_p << "\n"
                << "(dlnV/dlnp)_T   " << eq.state.dlnV_dlnp_T << "\n"
                << "sound speed     " << std::setprecision(2) << eq.state.a << " m/s\n\n"
                << "solver          " << eq.diagnostics.message << " after "
                << eq.diagnostics.restarts << " restart(s)\n"
                << std::scientific << std::setprecision(3)
                << "  element residual (relative) " << eq.diagnostics.element_residual_rel << "\n"
                << "  element residual (absolute) " << eq.diagnostics.element_residual_abs
                << " mol/kg\n"
                << "  Gibbs optimality residual   " << eq.diagnostics.gibbs_residual << "\n"
                << "  per-kg mass residual        " << eq.diagnostics.mass_residual << "\n"
                << "  state residual              " << eq.diagnostics.state_residual << "\n\n";
      std::cout << std::left << std::setw(10) << "species" << std::right << std::setw(16)
                << "mole frac" << std::setw(16) << "mass frac" << std::setw(18) << "mol/kg"
                << "\n";
      for (std::size_t k : order) {
        const auto i = static_cast<Eigen::Index>(k);
        if (X(i) < min_fraction) continue;
        std::cout << std::left << std::setw(10) << db[k].name() << std::right << std::scientific
                  << std::setprecision(6) << std::setw(16) << X(i) << std::setw(16) << Y(i)
                  << std::setw(18) << eq.state.n(i) << "\n";
      }
    }

    // --- export ----------------------------------------------------------
    Table t("equilibrium_composition");
    std::vector<std::string> names;
    std::vector<double> xs, ys, ns, mw;
    for (std::size_t k : order) {
      const auto i = static_cast<Eigen::Index>(k);
      names.push_back(db[k].name());
      xs.push_back(X(i));
      ys.push_back(Y(i));
      ns.push_back(eq.state.n(i));
      mw.push_back(db[k].molarMass());
    }
    t.addColumn("species", names);
    t.addColumn("mole_fraction", xs, "-");
    t.addColumn("mass_fraction", ys, "-");
    t.addColumn("moles_per_kg", ns, "mol/kg");
    t.addColumn("molar_mass", mw, "kg/mol");
    t.writeCsv(app::outputPath(cfg, "equilibrium.csv"));

    Json j = app::provenance(cfg);
    j["problem"] = Json(problem);
    j["mixture_ratio"] = Json(cfg.mixture_ratio);
    j["equivalence_ratio"] = Json(mix.equivalenceRatio());
    j["pressure"] = Json(cfg.chamber_pressure);
    j["reactant_enthalpy"] = Json(h0);
    j["temperature"] = Json(eq.state.T);
    j["molar_mass"] = Json(eq.state.M);
    j["density"] = Json(eq.state.rho);
    j["enthalpy"] = Json(eq.state.h);
    j["entropy"] = Json(eq.state.s);
    j["cp_frozen"] = Json(eq.state.cp_frozen);
    j["cp_equilibrium"] = Json(eq.state.cp_eff);
    j["cv_frozen"] = Json(eq.state.cv_frozen);
    j["cv_equilibrium"] = Json(eq.state.cv_eff);
    j["gamma_frozen"] = Json(eq.state.gamma_frozen);
    j["gamma_s"] = Json(eq.state.gamma_s);
    j["dlnV_dlnT_p"] = Json(eq.state.dlnV_dlnT_p);
    j["dlnV_dlnp_T"] = Json(eq.state.dlnV_dlnp_T);
    j["sound_speed"] = Json(eq.state.a);
    Json d = Json::object();
    d["iterations"] = Json(eq.diagnostics.iterations);
    d["restarts"] = Json(eq.diagnostics.restarts);
    d["element_residual_rel"] = Json(eq.diagnostics.element_residual_rel);
    d["element_residual_abs"] = Json(eq.diagnostics.element_residual_abs);
    d["gibbs_residual"] = Json(eq.diagnostics.gibbs_residual);
    d["mass_residual"] = Json(eq.diagnostics.mass_residual);
    d["state_residual"] = Json(eq.diagnostics.state_residual);
    d["convergence_history"] = Json::of(eq.diagnostics.history);
    j["diagnostics"] = d;
    j["composition"] = t.toJson();

    if (!flagPresent(cli, "no-transport")) {
      const TransportModel tr(db);
      const auto props = tr.mixture(X, eq.state.T, eq.state.cp_frozen);
      Json p = Json::object();
      p["viscosity"] = Json(props.viscosity);
      p["conductivity"] = Json(props.conductivity);
      p["prandtl"] = Json(props.prandtl);
      p["covered_mole_fraction"] = Json(props.covered_mole_fraction);
      j["transport"] = p;
      if (!cli.quiet)
        std::cout << "\ntransport (Chapman-Enskog, frozen cp)\n"
                  << std::scientific << std::setprecision(5)
                  << "  viscosity      " << props.viscosity << " Pa s\n"
                  << "  conductivity   " << props.conductivity << " W/(m K)\n"
                  << std::fixed << std::setprecision(5)
                  << "  Prandtl        " << props.prandtl << "\n"
                  << "  covered X      " << props.covered_mole_fraction << "\n";
    }
    // Optional mixture-ratio sweep of the whole composition.  The scalar
    // metric machinery in ignis_sweep cannot carry a composition vector, so the
    // per-species table is produced here.
    if (flagPresent(cli, "mr-sweep")) {
      const std::string arg = flagString(cli, "mr-sweep", "");
      double mr_lo = 0.0, mr_hi = 0.0;
      int mr_n = 0;
      if (std::sscanf(arg.c_str(), "%lf:%lf:%d", &mr_lo, &mr_hi, &mr_n) != 3 || mr_n < 2 ||
          !(mr_hi > mr_lo) || !(mr_lo > 0.0))
        throw ConfigError("--mr-sweep expects MIN:MAX:N with 0 < MIN < MAX and N >= 2, got '" +
                          arg + "'");
      Table sweep("composition_vs_mixture_ratio");
      std::vector<double> mrs, temps, molar, cstar_col;
      std::vector<std::vector<double>> frac(db.size());
      for (int i = 0; i < mr_n; ++i) {
        const double mr = mr_lo + (mr_hi - mr_lo) * i / (mr_n - 1);
        const PropellantMixture m2(ox, fu, mr, mix.oxidizerTemperature(), mix.fuelTemperature());
        const auto b2 = m2.elementMoles(db);
        const auto r2 = solver.hp(b2, m2.enthalpy(db), cfg.chamber_pressure);
        const NozzleFlow flow2(solver, CompositionModel::kEquilibrium,
                               ChamberReference{b2, r2.state});
        const auto X2 = r2.state.moleFractions();
        mrs.push_back(mr);
        temps.push_back(r2.state.T);
        molar.push_back(r2.state.M);
        cstar_col.push_back(flow2.cStarIdeal());
        for (std::size_t k = 0; k < db.size(); ++k)
          frac[k].push_back(X2(static_cast<Eigen::Index>(k)));
      }
      sweep.addColumn("mixture_ratio", mrs, "-");
      sweep.addColumn("temperature", temps, "K");
      sweep.addColumn("molar_mass", molar, "kg/mol");
      sweep.addColumn("c_star_ideal", cstar_col, "m/s");
      for (std::size_t k = 0; k < db.size(); ++k)
        sweep.addColumn("X_" + db[k].name(), frac[k], "-");
      sweep.writeCsv(app::outputPath(cfg, "composition_vs_mr.csv"));
      j["composition_vs_mixture_ratio"] = sweep.toJson();
      if (!cli.quiet)
        std::cout << "\nwrote " << app::outputPath(cfg, "composition_vs_mr.csv") << " ("
                  << mr_n << " mixture ratios)\n";
    }

    j.writeFile(app::outputPath(cfg, "equilibrium.json"));
    if (!cli.quiet)
      std::cout << "\nwrote " << app::outputPath(cfg, "equilibrium.csv") << "\n      "
                << app::outputPath(cfg, "equilibrium.json") << "\n";
    return 0;
  });
}
