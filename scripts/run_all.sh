#!/usr/bin/env bash
# Reproduce everything Ignis documents, from a clean checkout.
#
#   ./scripts/run_all.sh                 full reproduction
#   ./scripts/run_all.sh --quick         smaller Monte Carlo and optimisation
#   ./scripts/run_all.sh --stage build   run a single stage
#
# Stages: build tests methane hydrogen nozzle sweeps transient optimize mc
#         bench figures
#
# Every stage is independent; each prints what it produced.  The script exits
# non-zero on the first failure.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BUILD_DIR="${BUILD_DIR:-build}"
BIN="$BUILD_DIR/bin"
RESULTS="${RESULTS:-results}"
JOBS="${JOBS:-$( (command -v nproc >/dev/null && nproc) || echo 4)}"
QUICK=0
STAGES="build tests methane hydrogen nozzle sweeps transient optimize mc bench figures"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick) QUICK=1; shift ;;
    --stage) STAGES="$2"; shift 2 ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown option '$1' (try --help)" >&2; exit 2 ;;
  esac
done

MC_SAMPLES=2000
OPT_EVALS=2500
if [[ $QUICK -eq 1 ]]; then MC_SAMPLES=300; OPT_EVALS=600; fi

# The applications create their own output subdirectories, but the report files
# are written by `tee`, which opens its target as the pipeline starts.  Without
# this the first stage fails on a clean checkout where results/ does not exist.
mkdir -p "$RESULTS"

START_SECONDS=$SECONDS
banner() { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
has()    { [[ " $STAGES " == *" $1 "* ]]; }

# ---------------------------------------------------------------------------
if has build; then
  banner "configure and build  (jobs: $JOBS)"
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$BUILD_DIR" -j "$JOBS"
  "$BIN/ignis_engine" --version
fi

if has tests; then
  banner "test suite"
  ( cd "$BUILD_DIR" && ctest --output-on-failure -j "$JOBS" )
fi

if has methane; then
  banner "nominal LOX/methane engine"
  "$BIN/ignis_engine"      --config configs/methane_nominal.yaml --contour \
                           | tee "$RESULTS/methane_nominal_report.txt"
  "$BIN/ignis_equilibrium" --config configs/methane_nominal.yaml
  banner "LOX/methane, sea-level and vacuum operation of the vacuum nozzle"
  "$BIN/ignis_engine"      --config configs/methane_vacuum.yaml \
                           | tee "$RESULTS/methane_vacuum_report.txt"
  "$BIN/ignis_engine"      --config configs/methane_sea_level.yaml \
                           | tee "$RESULTS/methane_sea_level_report.txt"
fi

if has hydrogen; then
  banner "nominal LOX/hydrogen engine"
  "$BIN/ignis_engine"      --config configs/hydrogen_nominal.yaml \
                           | tee "$RESULTS/hydrogen_nominal_report.txt"
  "$BIN/ignis_equilibrium" --config configs/hydrogen_nominal.yaml --mr-sweep 3.0:8.0:26
fi

if has nozzle; then
  banner "altitude performance"
  "$BIN/ignis_nozzle" --config configs/methane_vacuum.yaml --points 161 --vacuum --contour \
                      | tee "$RESULTS/methane_vacuum_altitude_report.txt"
  "$BIN/ignis_nozzle" --config configs/methane_nominal.yaml --points 161 --vacuum
  banner "frozen against shifting equilibrium"
  "$BIN/ignis_nozzle" --config configs/methane_vacuum.yaml --points 41 --frozen \
                      --prefix m1vac_frozen --quiet
fi

if has sweeps; then
  banner "mixture-ratio sweep"
  "$BIN/ignis_sweep"       --config configs/mixture_ratio_sweep.yaml --mode sweep
  "$BIN/ignis_equilibrium" --config configs/mixture_ratio_sweep.yaml --mr-sweep 2.0:5.0:61 --quiet
  banner "expansion-ratio sweep"
  "$BIN/ignis_sweep"       --config configs/expansion_ratio_sweep.yaml --mode sweep
  banner "chamber-pressure sweep"
  "$BIN/ignis_sweep"       --config configs/chamber_pressure_sweep.yaml --mode sweep
  banner "altitude sweep"
  "$BIN/ignis_sweep"       --config configs/altitude_sweep.yaml --mode sweep
  banner "cooling-channel design space"
  "$BIN/ignis_sweep"       --config configs/cooling_sweep.yaml --mode sweep
fi

if has transient; then
  banner "startup and shutdown transient"
  "$BIN/ignis_transient" --config configs/startup_transient.yaml \
                         | tee "$RESULTS/startup_transient_report.txt"
  banner "integrator comparison (fixed-step RK4)"
  "$BIN/ignis_transient" --config configs/startup_transient.yaml --integrator rk4 --dt 2e-6 \
                         --prefix startup_rk4 --no-interp-check --quiet
fi

if has optimize; then
  banner "constrained design optimisation"
  "$BIN/ignis_sweep" --config configs/optimization.yaml --mode optimize \
                     | tee "$RESULTS/optimization_report.txt"
fi

if has mc; then
  banner "Monte Carlo campaign ($MC_SAMPLES samples)"
  "$BIN/ignis_mc" --config configs/monte_carlo.yaml --samples "$MC_SAMPLES" \
                  | tee "$RESULTS/monte_carlo_report.txt"
  banner "determinism check: one thread against many"
  "$BIN/ignis_mc" --config configs/monte_carlo.yaml --samples 200 --threads 1 \
                  --prefix mc_t1 --quiet
  "$BIN/ignis_mc" --config configs/monte_carlo.yaml --samples 200 --threads "$JOBS" \
                  --prefix mc_tN --quiet
  if cmp -s "$RESULTS/monte_carlo/mc_t1_mc_samples.csv" \
            "$RESULTS/monte_carlo/mc_tN_mc_samples.csv"; then
    echo "  single-thread and ${JOBS}-thread sample tables are byte-identical"
    # The two tables were the check; keeping identical copies of them is not
    # evidence of anything, so they go once the comparison has passed.
    rm -f "$RESULTS"/monte_carlo/mc_t1_* "$RESULTS"/monte_carlo/mc_tN_*
  else
    echo "  ERROR: Monte Carlo results depend on the thread count" >&2
    exit 1
  fi
fi

if has bench; then
  banner "benchmarks"
  "$BIN/ignis_bench" --output "$RESULTS/benchmarks" | tee "$RESULTS/benchmark_report.txt"
fi

if has figures; then
  banner "figures and animation"
  python3 python/make_figures.py --results "$RESULTS" --figures "$RESULTS/figures"
fi

banner "done"
ELAPSED=$((SECONDS - START_SECONDS))
printf 'results in %s/\n' "$RESULTS"
printf 'end-to-end wall time: %d min %02d s (%d s) on %d job(s)\n' \
       $((ELAPSED / 60)) $((ELAPSED % 60)) "$ELAPSED" "$JOBS"
printf 'stages run: %s\n' "$STAGES"
