# Architecture

Ignis is a layered C++17 library with six thin command-line front ends and a
Python post-processing package. Nothing in the library knows about files,
command lines or plots; nothing in the applications knows about
thermodynamics. Every module depends only on the modules below it, so the
dependency graph is a DAG and each layer can be tested on its own.

## Layer diagram

```
                        ┌──────────────────────────────────────────────┐
  applications          │ ignis_equilibrium  ignis_engine  ignis_nozzle│
  (apps/)               │ ignis_transient    ignis_sweep   ignis_mc    │
                        └───────────────────┬──────────────────────────┘
                                            │
                        ┌───────────────────┴──────────────────────────┐
  studies               │ optimize/  (Parameters, Sweep, Optimizer)    │
  (studies layer)       │ uncertainty/ (MonteCarlo, sensitivity)       │
                        └───────────────────┬──────────────────────────┘
                                            │
                        ┌───────────────────┴──────────────────────────┐
  assembly              │ engine/SteadyEngine  — one config → one      │
                        │ complete steady-state analysis               │
                        └───┬────────┬───────────┬──────────┬──────────┘
                            │        │           │          │
        ┌───────────────────┘        │           │          └──────────────┐
        │                            │           │                         │
  ┌─────┴──────┐        ┌────────────┴───┐  ┌────┴─────────┐   ┌──────────┴────┐
  │ combustion │        │ nozzle         │  │ thermal      │   │ cycle         │
  │ Propellant │        │ NozzleGeometry │  │ HeatTransfer │   │ FeedSystem    │
  │ Chamber    │        │ NozzleFlow     │  │ Regenerative │   │               │
  │            │        │ Atmosphere     │  │ CoolantFluid │   │               │
  └─────┬──────┘        └────────┬───────┘  └────┬─────────┘   └───────────────┘
        │                        │               │
        └──────────┬─────────────┴───────────────┘
                   │                                    ┌───────────────────┐
        ┌──────────┴───────────┐                        │ transient         │
        │ equilibrium          │◄───────────────────────┤ EquilibriumTable  │
        │ EquilibriumSolver    │                        │ TransientChamber  │
        └──────────┬───────────┘                        └───────────────────┘
                   │
        ┌──────────┴────────────────────────────────┐
        │ thermo                                    │
        │ Species  SpeciesDatabase  GasMixture      │
        │ GasState  Transport                       │
        └──────────┬────────────────────────────────┘
                   │
        ┌──────────┴────────────────────────────────┐
        │ core (Constants, Exceptions)   io (Config,│
        │ Json, Table, Cli)   Eigen   yaml-cpp      │
        └───────────────────────────────────────────┘
```

`io/` is a utility layer used by the applications and by the studies layer for
serialisation; the physics modules never touch it.

## Modules

| Module | Lines | Responsibility |
|---|---:|---|
| `core/` | 154 | Physical constants, unit and sign conventions, the exception hierarchy. Header-only. |
| `thermo/` | 1095 | NASA 7-coefficient species polynomials, the species database and its element matrix, per-kilogram gas mixtures with frozen-composition state functions, Chapman–Enskog transport. |
| `equilibrium/` | 813 | Constrained Gibbs-energy minimisation: (T,p), (h,p), (s,p) and (u,v) problems, equilibrium derivative properties (c<sub>p,eff</sub>, γ<sub>s</sub>, sound speed), element and energy residuals. |
| `combustion/` | 493 | Propellant definitions with liquid-state reference enthalpies (including RP-1 as CH<sub>1.95</sub>), mixture assembly by mixture ratio, chamber stagnation state, characteristic velocity. |
| `nozzle/` | 2142 | Rao-type bell contour generation, quasi-1D frozen and shifting-equilibrium marching, the throat solve, normal shocks, separation criteria, boundary-layer and kinetic corrections to thrust, US Standard Atmosphere 1976. |
| `thermal/` | 2786 | The integral turbulent boundary layer and Bartz's correlation for the hot-gas side, boundary-layer losses, Hatch & Papell film cooling, recovery temperature, wall conduction and radiation, tabulated reference-EOS coolant properties, the coupled regenerative-cooling march with tapered channels. |
| `kinetics/` | 915 | A gas-phase reaction mechanism (GRI-Mech 3.0's H/O/CO subset) with reverse rates by detailed balance, and the finite-rate march through the supersonic nozzle with its L-stable linearly implicit integrator. |
| `transient/` | 1180 | Pre-tabulated equilibrium EOS with tensor-product cubic-Hermite interpolation, 0-D chamber mass/energy ODEs, RK4 and adaptive Cash–Karp integration, conservation integrals. |
| `cycle/` | 1456 | Pressure-fed supply (injector, lines, tank pressures); real-fluid pumps, gas and expander turbines, generator and preburner combustion; gas-generator, staged-combustion and expander power balances. |
| `engine/` | 916 | `SteadyEngine` — the single entry point that turns one configuration into a complete steady-state analysis (chamber, nozzle, kinetics, boundary layer, performance, cooling, film, turbopump cycle, feed system). |
| `optimize/` | 1069 | Named parameter accessors over the configuration, N-dimensional sweeps, augmented-Lagrangian + Nelder–Mead constrained optimisation with Latin-hypercube multi-start. |
| `uncertainty/` | 639 | Distributions, deterministic per-sample seeding, threaded Monte Carlo, standardised regression coefficients, Spearman rank correlation, local elasticities. |
| `io/` | 1313 | YAML/JSON configuration loading with path-qualified error messages, JSON writing, CSV tables, the shared CLI parser. |
| `apps/` | 986 | The six executables. Each is argument parsing plus a call into the library plus output. |
| `tests/` | 6299 | 125 Catch2 test cases: unit, verification, validation against external references, integration. |

## The central abstraction

Everything above `thermo/` is expressed in terms of two objects:

* `GasMixture` — a species subset plus a composition vector `n` in **mol/kg of
  mixture**, so `Σ nⱼ Mⱼ = 1`. This is the CEA per-kilogram basis; it makes the
  element-conservation constraints linear with constant right-hand sides and
  makes every extensive property a specific property automatically.
* `EquilibriumResult` — composition, temperature, pressure and the full
  equilibrium derivative set, plus the residuals that prove the solve
  converged.

`NozzleFlow`, `TransientChamber` and `EquilibriumTable` all consume the same
`EquilibriumSolver`; there is exactly one implementation of the chemistry in
the repository.

## Data flow of a steady-state analysis

```
config.yaml
   │  io/Config
   ▼
EngineConfig ──► Propellant mixture (combustion/)
   │                     │
   │                     ▼
   │             Chamber stagnation state  ── equilibrium (h,p)
   │                     │                        │
   │                     ▼                        ▼
   │             c*_ideal, ṁ, throat solve   frozen/shifting properties
   │                     │
   │                     ▼
   ├──► NozzleGeometry ─► NozzleFlow march  ──► station table (p,T,ρ,u,M,A/At)
   │                             │
   │                             ├──► performance: F, Isp, CF, separation
   │                             │
   │                             ├──► finite-rate march (kinetics/) ──► kinetic efficiency
   │                             │
   │                             ▼
   │                     integral boundary layer, film  (thermal/)
   │                             │
   │                             ▼
   └──────────────────► regenerative jacket march (thermal/)
                                 │    ▲ jacket inlet pressure and temperature
                                 ▼    │ (the fuel pump's discharge)
                         coolant ΔT, Δp, T_wg, T_wc
                                 │
                                 ├──► boundary-layer losses ──► C_d, displaced exit, deficit
                                 │
                                 ▼
                         turbopump cycle (cycle/) ──► pump discharges, turbine,
                                 │                     generator flow, delivered Isp
                                 ▼
                         feed system (cycle/) ──► tank pressures
                                 │
                                 ▼
                         EngineResult ──► JSON + CSV + human-readable report
```

A single `EngineResult` carries every number the applications and the studies
layer report, so a sweep, an optimisation and a Monte Carlo campaign all run
the identical analysis path as a single nominal run — there is no reduced
"fast" model that could disagree with it.

## Error policy

Failures are exceptions, never sentinel values or silently clamped states:

| Exception | Raised when |
|---|---|
| `ConfigError` | A configuration key is missing, mistyped or out of range. Messages name the YAML path. |
| `RangeError` | A property is requested outside the validity range of its data (NASA polynomial limits, coolant table limits, material fit limits). |
| `ConvergenceError` | An iterative solve did not reach its tolerance, with the achieved residual in the message. |
| `InfeasibleError` | The physics has no solution for the requested inputs (coolant pressure reaches zero, the wall exceeds the melting point, a chamber mixture ratio outside the tabulated range). |

Sweeps, optimisation and Monte Carlo catch these per-point, count them by
class, and report them. A sweep line that reads `21 points, 2 failed` names
the first failures verbatim — a failed design is never written out as if it
had succeeded.

## Threading

Only three places use threads, all with deterministic results:

* `EquilibriumTable` build — independent grid points, order-independent.
* `MonteCarlo` — each sample seeds its own RNG with `SplitMix64(seed, index)`,
  so the sample set is byte-identical regardless of thread count. `run_all.sh`
  verifies this by diffing a 1-thread run against an N-thread run.
* `Sweep` — independent points, results written by index.

Everything else is single-threaded. The library holds no global mutable state.

## Build

CMake ≥ 3.16, C++17. `Eigen3` and `yaml-cpp` are found as system packages if
present and fetched with `FetchContent` otherwise; Catch2 v3 likewise, only
when `IGNIS_BUILD_TESTS=ON` (the default). Targets:

```
ignis                  static library (all modules)
ignis_equilibrium …    the six applications
ignis_bench            the benchmark harness
ignis_tests            the Catch2 suite, registered with CTest
```

`-DIGNIS_WERROR=ON` turns on `-Wall -Wextra -Wpedantic -Werror`; CI builds GCC
Release, GCC Debug and Clang Release with it.
