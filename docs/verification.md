# Verification

*Verification asks whether the equations are solved correctly. Validation —
whether the equations describe reality — is [`validation.md`](validation.md).*

Every number in this document was produced by running the suite. Reproduce it
with:

```bash
./scripts/build.sh
./scripts/test.sh          # the whole suite, through CTest
./scripts/validate.sh      # only the [verification] and [validation] cases,
                           # printing the measured errors
```

**Result: 75 test cases, 17,490 assertions, 0 failures** in 60 s (CTest, `-j4`), on GCC 13.3.0 and
Clang 18.1.3, Release and Debug, with `-Wall -Wextra -Wpedantic -Werror`.

---

## 1. Method-of-exact-solution checks

These compare the solver against a closed-form answer, so the only error is
the solver's own.

### 1.1 Quasi-1D nozzle against the analytic constant-γ nozzle

A synthetic perfect gas (constant c<sub>p</sub>, so γ = 1.4 exactly, built
through `SpeciesDatabase::fromSpecies`) is expanded through the real nozzle
solver and compared with the isentropic area–Mach relations.

| Quantity | Bound enforced | Largest observed |
|---|---|---|
| Sonic point p\*, T\*, c\* | 1e-8 rel | 4.5e-14 |
| p, T at M = 0.1, 0.3, 0.6, 0.9, 1.5, 2, 3, 4, 5 (both branches) | 2e-6 rel | 7.1e-12 |
| Recovered Mach from the area ratio | 2e-6 | 1.3e-10 |
| Mass flux ρuA − ṁ at ε = 1.2 … 40 | 1e-9 rel | below reporter precision |
| Stagnation enthalpy h + u²/2 − h₀ | 1e-12 rel | below reporter precision |
| γ<sub>frozen</sub> and γ<sub>s</sub> from the equilibrium derivative machinery | 1e-12 rel | 0 (exactly 1.4) |

The last row matters: the equilibrium-derivative code path is exercised on a
gas whose answer is known analytically, so a bug in ∂lnV/∂lnT or ∂lnV/∂lnp
would show up as γ<sub>s</sub> ≠ c<sub>p</sub>/c<sub>v</sub>.

### 1.2 Cylindrical wall conduction

The per-unit-gas-side-area wall resistance is checked against
r<sub>g</sub> ln((r<sub>g</sub>+t)/r<sub>g</sub>)/k — agreement to 1e-14
relative — and the resulting ΔT against q″ times that resistance to 1e-12.
The plane-wall limit is checked separately: at t/r<sub>g</sub> = 1e-3, 1e-4
and 1e-5 the cylindrical and plane forms differ by **less than
t/r<sub>g</sub>** in each case, which is what the series expansion requires
and is a stronger statement than a fixed tolerance.

### 1.3 Chamber filling with constant inflow

Driven with constant propellant flows and full heat release, the 0-D chamber
must settle on the state the steady relations give: ṁ<sub>out</sub> =
ṁ<sub>in</sub>, p<sub>c</sub> = ṁ c\*/A<sub>t</sub>, and T<sub>c</sub> equal to
the adiabatic flame temperature at the resulting state. All three hold to
1e-6, 1e-6 and 5e-3 relative respectively — the last is the interpolation
error of the tabulated EOS, measured in §2.4, not an integration error.
A separate case repeats the check at mid-burn on the shipped start-up
schedule (1e-4 on flow and pressure, 1e-3 on mixture ratio) and additionally
requires the chamber to be choked there.

### 1.3a The transient and the steady model as a cross-check

Those checks are internal to the transient. The stronger statement is against
the *other* solver. The steady model takes chamber pressure as an input and
returns mass flow; the transient takes mass flow as an input and returns
chamber pressure, through a completely separate path (tabulated equilibrium
EOS, ODE integration, orifice outflow). Running each on the other's answer:

| | |
|---|---:|
| `SteadyEngine` at p<sub>c</sub> = 5.500 MPa, O/F = 3.399267 | ṁ = 48.03612 kg/s |
| Transient commanded flows (37.12 + 10.92) | ṁ = 48.0400 kg/s |
| ⇒ steady model's implied plateau pressure | 5.50044 MPa |
| Transient's actual plateau pressure | **5.50212 MPa** |
| Disagreement | **0.031 %** |

This is *not* zero, and it should not be: the transient's ignition schedule
withholds a calibrated fraction of the heat of combustion
(`ignition.values` → 0.874), which is a different thing from the steady
model's η<sub>c\*</sub> = 0.96, and the tabulated EOS carries its own 0.13 %
interpolation error. 0.031 % is the residual after those two, and it is the
honest measure of how well two independent formulations of the same engine
agree.

### 1.4 Normal shock

The real-gas shock solver is run on the same synthetic perfect gas and
compared with the Rankine–Hugoniot relations for γ = 1.4 at upstream Mach 1.5,
2, 3 and 4. Pressure, density and temperature ratios and the downstream Mach
number match the closed-form relations to within 5e-6; mass, momentum and
total enthalpy across the jump are conserved to 1e-9. The entropy rise is
positive and the stagnation-pressure ratio below one at every Mach number, and
a subsonic upstream state raises `InfeasibleError` instead of returning a
spurious jump.

### 1.5 U.S. Standard Atmosphere 1976

| Altitude | Quantity | Ignis | Published (NOAA/NASA/USAF 1976) |
|---|---|---|---|
| 0 km | p | 101 325.0 Pa | 101 325 Pa |
| 0 km | T | 288.15 K | 288.15 K |
| 0 km | ρ | 1.224999 kg/m³ | 1.225 kg/m³ |
| 0 km | a | 340.2941 m/s | 340.294 m/s |
| 86 km | p | 0.3733805 Pa | 0.37338 Pa |
| 86 km | T | 186.946 K | 186.946 K |

Continuity across all seven layer joins is ≤ 3.2e-7 relative, and pressure
falls monotonically over the whole 0–86 km range.

---

## 2. Order-of-accuracy checks

### 2.1 Transient integrator

The fixed-step RK4 integrator is refined on a *smooth* chamber problem (fully
choked throughout, so the sub-critical orifice branch — which has a √
singularity at p<sub>c</sub> = p<sub>amb</sub> — is never crossed):

```
observed convergence orders: 4.14758, 4.05051, 4.17183
```

for mass, energy and pressure respectively — fourth order, as the method
requires.

> This test is written the way it is on purpose. Measured across the ignition
> transient of `configs/startup_transient.yaml`, the observed order is ≈ 1.5,
> because the chamber crosses p<sub>c</sub> = p<sub>amb</sub> once during
> start-up and the outflow relation is only C⁰ there. That is a property of
> the physical problem, not of the integrator, and
> [`limitations.md`](limitations.md) says so.

### 2.2 Adaptive against fixed-step

Cash–Karp RK4(5) at rtol 1e-9 and fixed-step RK4 at dt = 2e-6 s, on the full
start-up/shutdown problem:

```
largest pressure difference between integrators 2.43476e-07  (relative)
```

### 2.3 Grid convergence of the cooling march

The jacket is solved at 60, 120, 240, 480 and 960 segments. The change in
total heat load and in peak wall temperature falls monotonically with
refinement, and the final refinement (480 → 960) changes both by less than
0.2 %. The shipped 240 segments therefore sit inside the converged range.

### 2.4 Equilibrium property table

The transient's pre-tabulated EOS is checked at 200 random *interior* points
against a direct equilibrium solve at the same state, on the production grid
(41 × 97 × 13 = 51 701 nodes):

```
temperature : max 1.296e-01 %, rms 1.222e-02 %
pressure    : max 1.296e-01 %, rms 1.217e-02 %
molar mass  : max 7.120e-03 %
c*          : max 1.191e-02 %
```

Tensor-product cubic Hermite interpolation converges at fourth order in each
direction; linear interpolation on the same grid gives 4.8–6.2 %. Every
transient run prints this measurement, so the discretisation error is never
hidden.

---

## 3. Conservation and balance residuals

These are computed *after* each solve from the converged state, not assumed.

| Residual | Definition | Ignis-M1 | Ignis-H1 | Optimised design |
|---|---|---|---|---|
| Element balance | max<sub>i</sub> \|Σ<sub>j</sub> a<sub>ij</sub>n<sub>j</sub> − b<sub>i</sub>\| / b<sub>i</sub> | 1.14e-14 | 3.17e-15 | 9.09e-15 |
| Gibbs stationarity | max<sub>j</sub> \|μ<sub>j</sub>/RT − Σ<sub>i</sub> π<sub>i</sub>a<sub>ij</sub>\| | 1.42e-14 | 7.11e-15 | 1.78e-14 |
| Specified state (h) | \|h(n,T) − h<sub>0</sub>\| / \|h<sub>0</sub>\| | 7.36e-15 | 7.45e-15 | 6.53e-15 |
| Mass closure | Σ<sub>j</sub> n<sub>j</sub>M<sub>j</sub> − 1 | 1.11e-15 | −4.00e-15 | 6.22e-15 |
| Nozzle mass flux | max over stations of \|ρuA − ṁ\|/ṁ | 2.13e-14 | 4.37e-14 | 4.26e-14 |
| Nozzle stagnation enthalpy | max over stations of \|h + u²/2 − h<sub>0</sub>\|/\|h<sub>0</sub>\| | 2.94e-16 | 4.52e-16 | 1.52e-16 |
| Cooling global energy balance | \|ΣQ<sub>k</sub> − ṁ<sub>c</sub>Δh<sub>c</sub>\| / ΣQ<sub>k</sub> | 1.75e-10 | 4.83e-09 | 1.21e-09 |
| Cooling local flux consistency | max<sub>k</sub> spread between the gas-side, wall and coolant-side fluxes | 4.36e-10 | 3.02e-10 | 3.42e-10 |

For the start-up transient (`configs/startup_transient.yaml`, 802 output
samples over 0.8 s, Cash–Karp RK4(5) at rtol 1e-9):

| Residual | Definition | Measured |
|---|---|---|
| Mass conservation | \|∫(ṁ<sub>in</sub> − ṁ<sub>out</sub>)dt − Δm\| / ∫ṁ<sub>in</sub>dt | 5.708e-15 |
| Energy conservation | the same integral for the energy equation | 2.597e-13 |

These are *reported*, not assumed: the residuals are recomputed from the
converged state after every solve and printed by every application. A run
whose residuals are large is visible as such.

The element, Gibbs and state residuals are exposed as sweepable metrics
(`chamber.element_residual`, `chamber.gibbs_residual`,
`chamber.enthalpy_residual`), so they can be — and are — checked across an
entire parameter sweep rather than at a single point.

---

## 4. Consistency and invariance checks

| Check | What it enforces |
|---|---|
| Equilibrium from 12 random initial guesses (each species independently uniform on 1e-6 … 40 mol/kg) plus a default start | The solution is a property of the problem, not of the start point: T and M agree to 1e-8 relative and every mole fraction to 1e-8 absolute. |
| (s,p) and (T,v) problems started from a converged (h,p) state | T to 1e-9, s to 1e-10, p to 1e-9. |
| (u,v) problem started from the same state | T and p to 1e-7 — looser because the (u,v) problem is the least well conditioned of the four. |
| Equilibrium derivative terms against central finite differences of the solver itself | ∂lnV/∂lnT and ∂lnV/∂lnp to 2e-5, c<sub>p,eff</sub> to 5e-5, γ<sub>s</sub> to 1e-4 — these are the finite-difference floors, not solver error. |
| c<sub>p</sub> vs dh/dT and c<sub>p</sub>/T vs ds/dT for all 40 species over 200–6000 K | 1e-8 relative against a central difference. |
| Polynomial evaluation against an independent implementation of the same coefficients, 441 species/temperature pairs | Worst measured: c<sub>p</sub> 4.56e-12, h 4.55e-12, s 3.42e-12. |
| g = h − Ts for all 40 species | Exact by construction, checked to 1e-12. |
| NASA polynomial branch join at T<sub>mid</sub> | c<sub>p</sub> step < 2e-4, h and s steps < 2e-5 relative — that is the fitting precision of the source data, and it is measured rather than assumed. |
| Molar masses recomputed from the element matrix and atomic weights | 1e-6 relative for all 40 species. |
| A sweep point vs the equivalent single-point run | The sweep's O/F = 3.4 row reproduces a direct `SteadyEngine` run to 1e-12, and every metric column of a one-point sweep compares **bit-identical** (`==`, not a tolerance) to the matching row of the five-point sweep. The study layer runs the same code path; there is no reduced model that could disagree. |
| Monte Carlo at 1, 4 and 7 threads | Byte-identical sample matrices, NaN placement included. `scripts/run_all.sh` re-checks 1 vs 4 threads by `diff` on the written CSV for every full run. |
| Same seed twice / different seeds | Identical / different, as expected. |
| Drawn samples vs requested distributions | Over 4 000 samples: the normal input reproduces its mean to 1 %, its standard deviation to 6 % and has \|skewness\| < 0.15; truncated inputs never leave their bounds; the triangular input reproduces its analytic mean (low+mode+high)/3 to 1 % and stays inside its support; every output's order statistics satisfy min ≤ p5 ≤ p50 ≤ p95 ≤ p99 ≤ max. |
| Sensitivity ranking on a synthetic function with known analytic derivatives | Recovers the correct ordering and the correct signs of the standardised regression coefficients. |
| Optimiser result re-analysed from scratch | The reported objective and all constraint values reproduce, and the returned design is feasible. |

## 5. Error-path checks

A model that silently produces a plausible-looking number when it has failed
is worse than one that stops. These cases assert that Ignis stops, and that
the message is specific:

| Case | Expected |
|---|---|
| NASA polynomial evaluated outside 200–6000 K | `RangeError`, naming the species and the limit |
| Coolant property outside the table, or below the melting line | `RangeError`, naming the fluid, the state and the bound |
| Transport data missing for a species | `RangeError` — never a guessed Lennard-Jones pair |
| Equilibrium with an element that no species contains | `InfeasibleError` |
| Equilibrium driven to non-convergence (iteration cap) | `ConvergenceError` carrying the achieved residual |
| Nozzle expanded past the property range | `RangeError`, not a clamped state |
| Coolant pressure reaching zero in the jacket | `InfeasibleError` naming the axial position and the flow it cannot pass |
| Wall temperature above the material limit | `InfeasibleError` naming the material and the temperature |
| Transient mixture ratio outside the tabulated range | throws — it is **not** clamped, because clamping would make (u, MR) inconsistent |
| Malformed / missing / out-of-range configuration keys | `ConfigError` naming the YAML path and file |
| Unknown parameter or metric name in a study | `ConfigError` listing the valid names |
| Monte Carlo sample that raises | counted and classified by exception type, reported, never dropped |
| Invalid nozzle geometry (ε ≤ 1, negative radius, θ<sub>e</sub> > θ<sub>i</sub>, …) | `ConfigError` explaining which constraint failed |

The sweeps in `scripts/run_all.sh` exercise this for real: the
chamber-pressure sweep reports `21 points, 2 failed` and prints the two
`InfeasibleError` messages verbatim, and the cooling design-space sweep
reports `169 points, 52 failed`. Those failures are genuine — at 11 MPa the
shipped channel geometry cannot pass the flow — and they appear as failures,
not as optimistic numbers.

---

## 6. Reproducibility of the committed results

The `results/` tree in this repository was produced by `scripts/run_all.sh`.
To check that it can be reproduced, the branch was cloned into a fresh
directory, built from scratch and run end to end:

```bash
git clone --branch <branch> <url> /tmp/clean && cd /tmp/clean
rm -rf results && ./scripts/run_all.sh
```

Result: **exit 0 in 13 min 49 s**, 75/75 tests passed, all 16 figures written.
Comparing the fresh output against the committed one, with the two provenance
lines (git hash, data path) excluded:

| Compared | Files | Files with any difference | What the differences were |
|---|---:|---:|---|
| Human-readable reports | 9 | 4 | one `wall time …` line each |
| CSV tables | 33 | 4 | the four benchmark timing tables |

Every physics number — composition, temperature, c\*, thrust, I<sub>sp</sub>,
heat flux, wall temperature, pressure drop, every sweep point, every Monte
Carlo sample, every residual — is **bit-identical**. The only things that move
between runs are measured times, which is what should move. This was checked
twice, at two different commits.

This check also found a real bug: `run_all.sh` let the applications create
their own output subdirectories but wrote the reports through `tee`, which
opens its target as the pipeline starts. On a machine where `results/` already
existed it worked; from a clean checkout the first stage died with
`tee: results/methane_nominal_report.txt: No such file or directory`. That is
the entire argument for running your own instructions against a fresh clone.

---

## 7. Test inventory

| File | Cases | Focus |
|---|---:|---|
| `tests/unit/test_thermo.cpp` | 10 | NASA polynomials, database, mixtures, range errors |
| `tests/unit/test_equilibrium.cpp` | 8 | Conservation, stationarity, guess independence, trends, derivatives, error paths |
| `tests/unit/test_nozzle.cpp` | 6 | Contour validity, analytic nozzle, thrust decomposition, regime classification |
| `tests/unit/test_transport_atmosphere.cpp` | 5 | Collision integrals, Chapman–Enskog, Wilke, USSA-1976 |
| `tests/unit/test_thermal.cpp` | 8 | Conduction, recovery temperature, Bartz scaling, radiation, coolant tables |
| `tests/unit/test_cooling.cpp` | 6 | Energy balance, monotonicity, grid convergence, design trends, failures |
| `tests/unit/test_transient.cpp` | 7 | Table accuracy, steady-state agreement, analytic filling, order, adaptivity, shutdown |
| `tests/unit/test_io.cpp` | 7 | Every shipped config, defaults, error messages, parameter registry, JSON/CSV |
| `tests/unit/test_uncertainty.cpp` | 7 | Thread invariance, seeding, distributions, sensitivity, failure accounting |
| `tests/validation/test_validation.cpp` | 5 | NASA CEA and Cantera comparisons — see [`validation.md`](validation.md) |
| `tests/integration/test_engine.cpp` | 6 | Both engines end to end, altitude trends, sweep/single-point identity, optimisation feasibility, feed system |
| **Total** | **75** | **17 490 assertions, 60 s at `-j4`** |
