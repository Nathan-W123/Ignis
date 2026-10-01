# Benchmarks

Every number here was produced by `ignis_bench` on the machine described
below. Nothing is estimated, extrapolated or copied from another machine.

```bash
./scripts/build.sh
./build/bin/ignis_bench -o results/benchmarks       # ~4 minutes
```

The harness writes `results/benchmarks/benchmarks.{csv,json}` alongside the
console report, and `scripts/run_all.sh` runs it as the `bench` stage.

## Test configuration

| | |
|---|---|
| CPU | Intel(R) Xeon(R) Processor @ 2.80 GHz |
| Hardware threads | 4 |
| Compiler | GCC 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) |
| Build type | `Release` (`-O3 -DNDEBUG`) |
| Eigen | 3.4.0 |
| yaml-cpp | 0.8.0 |
| Catch2 | 3.4.0 (system) |
| CMake | 3.28.3 |
| Python | NumPy 2.4.6, pandas 3.0.5, Matplotlib 3.11.2 |
| Ignis | 1.0.0, `v1.0.0-28-gfa97611`, a clean tree |

Timings are the **mean** of *n* repeats after a warm-up; the min and max of the
same set are reported so the spread is visible. The benchmark harness
re-checks the physics residuals while it times, so a fast-but-wrong run cannot
pass unnoticed.

## Core solver

| Benchmark | Mean | Min | Max | n | Secondary |
|---|---:|---:|---:|---:|---|
| Adiabatic (h,p) equilibrium, 26 species | **0.0588 ms** | 0.0531 | 0.4237 | 500 | 24.0 mean Newton iterations |
| Isentropic (s,p) equilibrium, 26 species | **0.1433 ms** | 0.1219 | 0.8512 | 500 | — |

Worst residuals over those 500 adiabatic solves: element **2.21e-14**, Gibbs
**5.68e-14**, enthalpy **4.12e-14**. Speed is not being bought with tolerance.

### Scaling with species count

| Species | Mean time | Relative to 8 species |
|---:|---:|---:|
| 8 | 0.0115 ms | 1.00× |
| 12 | 0.0266 ms | 2.31× |
| 18 | 0.0311 ms | 2.70× |
| 26 | 0.0611 ms | 5.31× |

The Newton system is (E+1)×(E+1) at fixed temperature and (E+2)×(E+2) when the
temperature is unknown — 4×4 or 5×5 for a C/H/O system, **independent of the
species count**, because the per-species unknowns are eliminated analytically.
(These are adiabatic solves, so 5×5.) The growth above is the O(N·E) assembly
of that system and the per-species Gibbs evaluations, not a growing linear
solve. The measured exponent is 1.42 over this range, consistent with linear
assembly plus a fixed overhead.

## Nozzle

The inviscid shifting-equilibrium nozzle alone: no wall boundary layer, no
finite-rate march, no jacket and no cycle. Those are timed in the next section.

| Benchmark | Mean | Min | Max | n |
|---|---:|---:|---:|---:|
| Chamber + throat + exit state (no profile) | **1.50 ms** | 1.39 | 1.97 | 50 |

### Scaling with station count

| Stations | Mean time | ms per station |
|---:|---:|---:|
| 50 | 5.03 ms | 0.101 |
| 100 | 7.62 ms | 0.076 |
| 200 | 11.35 ms | 0.057 |
| 400 | 20.74 ms | 0.052 |
| 800 | 37.76 ms | 0.047 |
| 1600 | 68.05 ms | 0.043 |

Linear in station count, as a marching scheme must be, with the fixed chamber
and throat solve amortised away by a few hundred stations. Each station is an
area-ratio inversion: a Newton iteration on the exact slope of ρu along the
isentrope ([`theory.md` §7](theory.md)), three or four warm-started
shifting-equilibrium (s,p) solves. The marginal cost is about 0.04 ms a
station. Before that iteration replaced a bracketing secant search of about 35
solves, a station cost about 0.4 ms and the chamber + throat + exit state
7.3 ms.

## End-to-end analysis

`configs/methane_nominal.yaml` in full, then with one stage left out at a time:

| Benchmark | Mean | Min | Max | n |
|---|---:|---:|---:|---:|
| Everything: chamber, 400-station profile, boundary layer, finite-rate march, jacket with film, gas-generator cycle, feed | **351.3 ms** | 325.5 | 395.5 | 12 |
| Without the finite-rate march | 218.1 ms | 199.2 | 231.8 | 12 |
| Without the boundary-layer losses | 331.2 ms | 296.6 | 377.2 | 12 |
| Without the jacket, cycle and feed | 190.2 ms | 173.3 | 239.4 | 12 |

Worst residuals during those runs: cooling energy balance **7.49e-14**, local
flux **1.95e-09**, nozzle mass flow **1.29e-13**.

Read as differences: the finite-rate march costs about 133 ms, 38 % of a
complete analysis; the jacket with its coupled boundary-layer passes, the
cycle and the feed together about 161 ms; the boundary-layer losses about
20 ms. The chamber, nozzle and profile are the remaining few tens of
milliseconds. A study that does not read the profile turns it off
(`sample_profile = false`), and neither the shipped trade study nor the Monte
Carlo campaign runs the finite-rate march or the cycle.

## Transient

| Benchmark | Mean | n | Secondary |
|---|---:|---:|---|
| Equilibrium table build (41 × 97 × 13) | **1439.2 ms** | 1 | 51 701 equilibrium solves |
| Adaptive RK4(5), 0.8 s of physical time | **1452.3 ms** | 3 | 4624 accepted steps |
| Fixed-step RK4 at dt = 2 µs, 0.8 s | **72 495 ms** | 1 | 400 155 steps |

Conservation residuals on the adaptive run: mass **5.71e-15**, energy
**2.60e-13**. Table interpolation error over 100 random interior points:
max |ΔT|/T **1.30e-1 %**, rms **1.69e-2 %**.

Two things are worth reading off this table:

* The table build costs 1.44 s for 51 701 equilibrium solves — **27.8 µs per
  solve** across 4 threads, against 58.8 µs single-threaded, a 2.1× speed-up.
* The adaptive integrator is **50× faster** than fixed-step RK4 at the step
  size needed for the same accuracy, and the two agree to 2.4e-7 relative in
  chamber pressure ([`verification.md` §2.2](verification.md)). Tabulating the
  equilibrium EOS is what makes either possible. Arithmetic, not a measurement:
  the fixed-step run evaluates 4 right-hand sides per step, so 400 155 steps
  would need 1.6 M equilibrium solves; at the measured 58.8 µs each that is
  **94 s of chemistry alone**, more than the 72.5 s the whole tabulated run
  takes. The table is built once in 1.4 s and read 1.6 M times.

## Monte Carlo throughput

400 samples of the full `configs/monte_carlo.yaml` analysis (integral boundary
layer, film, jacket; no finite-rate march or cycle):

| Threads | Wall time | Samples/s | Speed-up | Efficiency |
|---:|---:|---:|---:|---:|
| 1 | 83.45 s | 4.79 | 1.00× | 100 % |
| 2 | 43.58 s | 9.18 | 1.91× | 96 % |
| 4 | 24.92 s | 16.05 | **3.35×** | 84 % |

On 4 hardware threads the production campaign (2000 samples) completes in
**102.2 s at 19.6 samples/s**, slightly faster per sample than the benchmark
because the benchmark includes its own set-up. Parallel efficiency falls to
84 % at 4 threads mainly because the samples are not equally expensive: a
sample whose jacket nearly fails takes several times longer than a nominal
one, so the tail straggles.

Determinism costs nothing here: each sample seeds its own RNG from
`SplitMix64(seed, index)`, so the 1-thread and 4-thread sample matrices are
byte-identical, which `run_all.sh` verifies with `cmp` on every full run.

## Study wall times

Measured during the regeneration of `results/` (`scripts/run_all.sh`):

| Study | Points / samples | Wall time |
|---|---:|---:|
| Mixture-ratio sweep | 31 | 1.64 s |
| Expansion-ratio sweep | 160 | 1.76 s |
| Chamber-pressure sweep | 21 | 1.10 s |
| Altitude sweep | 81 | 0.87 s |
| Cooling design space | 169 (31 infeasible) | 7.15 s |
| Start-up transient (RK4(5)) | 0.8 s physical | 1.43 s |
| Constrained ascent trade study | 2500 evaluations, 12 starts (64 failed analyses) | 418.5 s |
| Monte Carlo campaign | 2000 samples | 102.2 s |
| Benchmark suite | — | 244.6 s |
| Test suite (CTest, 4 jobs) | 125 cases | 95.9 s |

The expansion-ratio sweep is fast because only the nozzle is re-solved; the
chamber state is shared. The cooling sweep is slow per point because each
point runs a full jacket march with its coupled boundary-layer passes.

## Build and reproduction time

| Step | Measured |
|---|---:|
| Clean build, 60 translation units, `-j4`, Release | **56 s** |
| Full reproduction from a fresh clone, build included (`scripts/run_all.sh`) | **19 min 58 s** |

The configure step is fast because Eigen, yaml-cpp and Catch2 are present as
system packages on this machine. Where they are not, CMake fetches them with
`FetchContent`, and the configure step then costs whatever the three clones
cost on the network in question — typically a minute or two, once.

`scripts/run_all.sh` prints its own end-to-end wall time on completion, so the
figure above is re-measured on every run rather than being a claim in a
document.

`--quick` shortens the three stages that dominate a full run — 300 Monte Carlo
samples instead of 2000, 600 optimiser evaluations instead of 2500, a reduced
benchmark pass, and no fixed-step RK4 comparison. It is what CI uses.

### The test suite

| | |
|---|---:|
| Full suite, CTest at `-j4` | **95.9 s** (125 cases, 21 359 assertions) |
| Longest single case (`constrained optimisation respects its constraints`) | 94.2 s |
| Quick check, `ignis_tests "~[slow]"`, one process | 38.9 s (121 cases) |

The trade-study case runs 400 complete engine evaluations, so it sets the
suite's wall time at `-j4` and follows the engine's own cost. When the
boundary layer, its losses and the finite-rate march went in, an evaluation
became several times more expensive and the suite took 16.5 minutes in CI.
The Newton area-ratio search and two caches brought it back, moving the
Ignis-M1's thrust, impulses, kinetic efficiency and wall temperatures by no
more than 10⁻⁹ ([`verification.md` §7](verification.md)).

That case once took 257 s for a different reason: its six `SECTION`s made
Catch2 re-run the test body, and therefore the whole optimisation, once per
section. Two Monte Carlo cases had the same problem. The checks now share one
campaign each, with the same assertions.
