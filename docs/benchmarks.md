# Benchmarks

Every number here was produced by `ignis_bench` on the machine described
below. Nothing is estimated, extrapolated or copied from another machine.

```bash
./scripts/build.sh
./build/bin/ignis_bench -o results/benchmarks       # ~3 minutes
```

The harness writes `results/benchmarks/benchmarks.{csv,json}` alongside the
console report, and `scripts/run_all.sh` runs it as the `bench` stage.

## Test configuration

| | |
|---|---|
| CPU | Intel(R) Xeon(R) Processor @ 2.80 GHz |
| Hardware threads | 4 |
| Compiler | GCC 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) |
| Build type | `Release` (`-O3 -DNDEBUG`), `-DIGNIS_WERROR=ON` |
| Eigen | 3.4.0 |
| yaml-cpp | 0.8.0 |
| Catch2 | 3.4.0 (system) / 3.5.2 (fetched) |
| CMake | 3.28.3 |
| Python | NumPy 2.4.6, pandas 3.0.5, Matplotlib 3.11.2 |
| Ignis | 1.0.0 |

Timings are the **mean** of *n* repeats after a warm-up; the min and max of the
same set are reported so the spread is visible. The benchmark harness
re-checks the physics residuals while it times, so a fast-but-wrong run cannot
pass unnoticed.

## Core solver

| Benchmark | Mean | Min | Max | n | Secondary |
|---|---:|---:|---:|---:|---|
| Adiabatic (h,p) equilibrium, 26 species | **0.0577 ms** | 0.0533 | 0.1654 | 500 | 24.0 mean Newton iterations |
| Isentropic (s,p) equilibrium, 26 species | **0.1404 ms** | 0.1238 | 0.4558 | 500 | — |

Worst residuals over those 500 adiabatic solves: element **2.21e-14**, Gibbs
**5.68e-14**, enthalpy **4.12e-14**. Speed is not being bought with tolerance.

### Scaling with species count

| Species | Mean time | Relative to 8 species |
|---:|---:|---:|
| 8 | 0.0132 ms | 1.00× |
| 12 | 0.0241 ms | 1.83× |
| 18 | 0.0308 ms | 2.34× |
| 26 | 0.0586 ms | 4.44× |

The Newton system is (E+1)×(E+1) at fixed temperature and (E+2)×(E+2) when the
temperature is unknown — 4×4 or 5×5 for a C/H/O system, **independent of the
species count**, because the per-species unknowns are eliminated analytically.
(These are adiabatic solves, so 5×5.) The growth above is the O(N·E) assembly
of that system and the per-species Gibbs evaluations, not a growing linear
solve. The measured exponent is 1.27 over this range, consistent with linear
assembly plus a fixed overhead.

## Nozzle

| Benchmark | Mean | Min | Max | n |
|---|---:|---:|---:|---:|
| Chamber + throat + exit state (no profile, no cooling) | **7.99 ms** | 7.56 | 12.02 | 50 |

### Scaling with station count

| Stations | Mean time | ms per station |
|---:|---:|---:|
| 50 | 33.02 ms | 0.660 |
| 100 | 49.89 ms | 0.499 |
| 200 | 88.13 ms | 0.441 |
| 400 | 163.51 ms | 0.409 |
| 800 | 319.46 ms | 0.399 |
| 1600 | 625.99 ms | 0.391 |

Linear in station count, as a marching scheme must be, with the fixed chamber
and throat solve amortised away by 400 stations. Each station is a full
shifting-equilibrium (s,p) solve plus a root-find on the area ratio, which is
why a station costs ≈ 0.4 ms while a bare equilibrium solve costs 0.06 ms.

## End-to-end analysis

| Benchmark | Mean | Min | Max | n |
|---|---:|---:|---:|---:|
| Chamber + nozzle + 400-station profile + cooling + feed | **263.24 ms** | 258.55 | 275.88 | 12 |
| The same without the regenerative jacket | **161.71 ms** | 160.68 | 163.14 | 12 |

Worst residuals during those runs: cooling energy balance **1.75e-10**, local
flux **4.36e-10**, nozzle mass flow **1.68e-13**.

So the 240-segment regenerative jacket costs about 102 ms, or 39 % of a
complete analysis; the axial profile costs about 154 ms of the rest. A study
that does not read the profile turns it off (`sample_profile = false`), which
is why the Monte Carlo campaign runs far faster per sample than 263 ms would
suggest.

## Transient

| Benchmark | Mean | n | Secondary |
|---|---:|---:|---|
| Equilibrium table build (41 × 97 × 13) | **1215.8 ms** | 1 | 51 701 equilibrium solves |
| Adaptive RK4(5), 0.8 s of physical time | **1237.5 ms** | 3 | 4624 accepted steps |
| Fixed-step RK4 at dt = 2 µs, 0.8 s | **65 626 ms** | 1 | 400 155 steps |

Conservation residuals on the adaptive run: mass **5.71e-15**, energy
**2.60e-13**. Table interpolation error over 100 random interior points:
max |ΔT|/T **1.30e-1 %**, rms **1.69e-2 %**.

Two things are worth reading off this table:

* The table build costs 1.22 s for 51 701 equilibrium solves — **23.5 µs per
  solve** across 4 threads, against 57.7 µs single-threaded, a 2.45× speed-up.
* The adaptive integrator is **53× faster** than fixed-step RK4 at the step
  size needed for the same accuracy, and the two agree to 2.4e-7 relative in
  chamber pressure ([`verification.md` §2.2](verification.md)). Tabulating the
  equilibrium EOS is what makes either possible. Arithmetic, not a measurement:
  the fixed-step run evaluates 4 right-hand sides per step, so 400 155 steps
  would need 1.6 M equilibrium solves; at the measured 57.7 µs each that is
  **92 s of chemistry alone**, more than the 65.6 s the whole tabulated run
  takes. The table is built once for 1.2 s and read 1.6 M times.

## Monte Carlo throughput

400 samples of the full `configs/monte_carlo.yaml` analysis:

| Threads | Wall time | Samples/s | Speed-up | Efficiency |
|---:|---:|---:|---:|---:|
| 1 | 46.58 s | 8.59 | 1.00× | 100 % |
| 2 | 24.75 s | 16.16 | 1.88× | 94 % |
| 4 | 13.70 s | 29.20 | **3.40×** | 85 % |

On 4 hardware threads the production campaign (2000 samples) completes in
**58.20 s at 34.4 samples/s** — slightly faster per sample than the benchmark
because the benchmark includes its own set-up. Parallel efficiency falls to
85 % at 4 threads mainly because the samples are not equally expensive: a
sample whose jacket nearly fails takes several times longer than a nominal
one, so the tail straggles.

Determinism costs nothing here: each sample seeds its own RNG from
`SplitMix64(seed, index)`, so the 1-thread and 4-thread sample matrices are
byte-identical, which `run_all.sh` verifies with `cmp` on every full run.

## Study wall times

Measured during the full reproduction (`scripts/run_all.sh`):

| Study | Points / samples | Wall time |
|---|---:|---:|
| Mixture-ratio sweep | 31 | 0.91 s |
| Expansion-ratio sweep | 160 | 0.10 s |
| Chamber-pressure sweep | 21 (2 infeasible) | 0.63 s |
| Altitude sweep | 81 | 0.05 s |
| Cooling design space | 169 (52 infeasible) | 4.48 s |
| Start-up transient (RK4(5)) | 0.8 s physical | 1.30 s |
| Constrained optimisation | 2500 evaluations (116 failed) | 245.06 s |
| Monte Carlo campaign | 2000 samples | 58.20 s |
| Benchmark suite | — | 169.03 s |
| Test suite (CTest, 4 jobs) | 75 cases | 257.16 s |

The expansion-ratio sweep is fast because only the nozzle is re-solved; the
chamber state is shared. The cooling sweep is slow per point because each
point runs a full jacket march.

## Build and reproduction time

| Step | Measured |
|---|---:|
| Configure from an empty build directory | **0.3 s** |
| Clean build, 45 translation units, `-j4`, Release + `-Werror` | **36 s** |
| Full reproduction from a warm build (`scripts/run_all.sh`) | **16 min 54 s** |

The configure step is fast because Eigen, yaml-cpp and Catch2 are present as
system packages on this machine. Where they are not, CMake fetches them with
`FetchContent`, and the configure step then costs whatever the three clones
cost on the network in question — typically a minute or two, once.

`scripts/run_all.sh` prints its own end-to-end wall time on completion, so the
figure above is re-measured on every run rather than being a claim in a
document. `--quick` (300 Monte Carlo samples, 600 optimiser evaluations) cuts
it to a few minutes and is what CI uses.
