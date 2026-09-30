# Ignis documentation

| Document | Read it for |
|---|---|
| [`theory.md`](theory.md) | The complete mathematical formulation — 17 sections, every equation and every correlation with its source |
| [`architecture.md`](architecture.md) | Module layout, dependency graph, data flow, error policy, threading model |
| [`configuration.md`](configuration.md) | Every configuration key, and every named parameter and metric a study can address |
| [`verification.md`](verification.md) | Are the equations solved correctly? Exact solutions, order of accuracy, conservation residuals, error paths |
| [`validation.md`](validation.md) | Do the equations describe reality? NASA CEA, Cantera and reference-EOS comparisons, with measured errors |
| [`benchmarks.md`](benchmarks.md) | Measured timings, scaling with problem size, parallel efficiency, build and reproduction times |
| [`limitations.md`](limitations.md) | What the model cannot do, why, and what would have to change |
| [`explorer.md`](explorer.md) | The desktop Engine Explorer: layout, charts, constraint panel, colour contract |
| [`extending.md`](extending.md) | How to add a species, a propellant, a correlation, a metric, a figure |

Start with the [root README](../README.md) for what Ignis is and what it
produces; come here for why the numbers are what they are.

Every number quoted in `verification.md`, `validation.md` and `benchmarks.md`
was produced by running the code on the machine described in `benchmarks.md`.
Reproduce them with:

```bash
./scripts/build.sh
./scripts/test.sh        # the whole suite
./scripts/validate.sh    # the V&V cases, printing their measured errors
./build/bin/ignis_bench -o results/benchmarks
./scripts/run_all.sh     # everything, including the figures
```
